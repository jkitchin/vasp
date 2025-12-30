"""Tests for the main Vasp calculator class."""

import os
import pickle

import numpy as np
import pytest

from vasp import CalculationResult, Vasp
from vasp.exceptions import (
    VaspQueued,
    VaspRunning,
    VaspSubmitted,
)
from vasp.runners import JobState, LocalRunner, MockResults, MockRunner


class TestVaspInit:
    """Test Vasp calculator initialization."""

    def test_basic_init(self, temp_dir, si_atoms):
        """Test basic calculator creation."""
        calc = Vasp(label=temp_dir, atoms=si_atoms)

        assert calc.directory == os.path.abspath(temp_dir)
        assert calc.parameters['xc'] == 'PBE'
        assert calc.parameters['kpts'] == (1, 1, 1)
        assert isinstance(calc.runner, LocalRunner)

    def test_init_with_parameters(self, temp_dir, si_atoms):
        """Test calculator with custom parameters."""
        calc = Vasp(
            label=temp_dir,
            atoms=si_atoms,
            xc='HSE06',
            encut=500,
            kpts=(4, 4, 4),
            ismear=0,
            sigma=0.05,
        )

        assert calc.parameters['xc'] == 'HSE06'
        assert calc.parameters['encut'] == 500
        assert calc.parameters['kpts'] == (4, 4, 4)
        assert calc.parameters['ismear'] == 0
        assert calc.parameters['sigma'] == 0.05
        # HSE06 sets these automatically
        assert calc.parameters['lhfcalc'] is True
        assert calc.parameters['hfscreen'] == 0.2

    def test_init_with_runner(self, temp_dir, si_atoms, mock_runner):
        """Test calculator with custom runner."""
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=mock_runner)

        assert calc.runner is mock_runner

    def test_xc_defaults(self, temp_dir, si_atoms):
        """Test that XC functional sets appropriate defaults."""
        # Test PBE0
        calc = Vasp(label=temp_dir, atoms=si_atoms, xc='PBE0')
        assert calc.parameters['lhfcalc'] is True
        assert calc.parameters['aexx'] == 0.25

        # Test SCAN meta-GGA
        calc = Vasp(label=temp_dir, atoms=si_atoms, xc='SCAN')
        assert calc.parameters['metagga'] == 'SCAN'

        # Test vdW-DF
        calc = Vasp(label=temp_dir, atoms=si_atoms, xc='optPBE-vdW')
        assert calc.parameters['luse_vdw'] is True
        assert calc.parameters['gga'] == 'OR'


class TestAtomSorting:
    """Test atom sorting for VASP ordering."""

    def test_sorting_same_species(self, temp_dir, si_atoms):
        """Test sorting when all atoms are same species."""
        calc = Vasp(label=temp_dir, atoms=si_atoms)

        assert calc.sort == [0, 1]
        assert calc.resort == [0, 1]

    def test_sorting_multiple_species(self, temp_dir):
        """Test sorting with multiple species."""
        from ase import Atoms

        # Create mixed species structure
        atoms = Atoms(
            'SiCSiC',
            positions=[
                [0, 0, 0],
                [1, 0, 0],
                [2, 0, 0],
                [3, 0, 0],
            ],
            cell=[10, 10, 10],
            pbc=True
        )

        calc = Vasp(label=temp_dir, atoms=atoms)

        # Should group Si together, then C together
        assert calc.sort == [0, 2, 1, 3]  # Si at 0,2 then C at 1,3
        assert calc.resort == [0, 2, 1, 3]  # Maps sorted back to original

    def test_resort_inverse(self, temp_dir):
        """Test that resort is the inverse of sort."""
        from ase import Atoms

        atoms = Atoms(
            'FeCoNiFeCo',
            positions=np.random.rand(5, 3) * 10,
            cell=[10, 10, 10],
            pbc=True
        )

        calc = Vasp(label=temp_dir, atoms=atoms)

        # Apply sort then resort should give identity
        original = list(range(len(atoms)))
        sorted_idx = [original[i] for i in calc.sort]
        resorted = [sorted_idx[i] for i in calc.resort]

        assert resorted == original


class TestCalculate:
    """Test calculation execution."""

    def test_calculate_complete(self, vasp_calc):
        """Test calculation that completes immediately."""
        vasp_calc.calculate()

        assert 'energy' in vasp_calc.results
        assert 'forces' in vasp_calc.results
        assert vasp_calc.results['energy'] == pytest.approx(-10.84274516)

    def test_get_potential_energy(self, vasp_calc):
        """Test getting potential energy."""
        energy = vasp_calc.get_potential_energy()
        assert energy == pytest.approx(-10.84274516)

    def test_get_forces(self, vasp_calc):
        """Test getting forces."""
        forces = vasp_calc.get_forces()

        assert forces.shape == (2, 3)
        assert forces[0, 0] == pytest.approx(0.012345)
        assert forces[1, 0] == pytest.approx(-0.012345)

    def test_get_stress(self, vasp_calc):
        """Test getting stress tensor."""
        stress = vasp_calc.get_stress()

        assert stress.shape == (6,)
        # Mock stress is [1.234, ...] in kBar, converted to eV/A^3 by * -0.1
        assert stress[0] == pytest.approx(-0.1234)

    def test_potential_energy_property(self, vasp_calc):
        """Test potential_energy property."""
        energy = vasp_calc.potential_energy
        assert energy == pytest.approx(-10.84274516)

    def test_forces_property(self, vasp_calc):
        """Test forces property."""
        forces = vasp_calc.forces
        assert forces.shape == (2, 3)


class TestAsyncWorkflow:
    """Test async/non-blocking workflow with exceptions."""

    def test_submitted_exception(self, temp_dir, si_atoms):
        """Test that VaspSubmitted is raised on first submit."""
        runner = MockRunner(
            results=MockResults(energy=-10.0),
            state_sequence=[JobState.SUBMITTED]
        )
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=runner)

        with pytest.raises(VaspSubmitted) as exc_info:
            calc.calculate()

        # Check that jobid has the expected format
        assert exc_info.value.jobid.startswith('mock-job-')

    def test_queued_exception(self, temp_dir, si_atoms):
        """Test that VaspQueued is raised when job is queued."""
        runner = MockRunner(
            results=MockResults(energy=-10.0),
            state_sequence=[JobState.QUEUED]
        )
        # Simulate already submitted
        runner._call_count = 0

        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=runner)

        with pytest.raises(VaspQueued):
            calc.calculate()

    def test_running_exception(self, temp_dir, si_atoms):
        """Test that VaspRunning is raised when job is running."""
        runner = MockRunner(
            results=MockResults(energy=-10.0),
            state_sequence=[JobState.RUNNING]
        )

        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=runner)

        with pytest.raises(VaspRunning):
            calc.calculate()

    def test_full_async_sequence(self, temp_dir, si_atoms, async_mock_runner):
        """Test full async workflow through multiple states."""
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=async_mock_runner)

        # First call - submitted
        with pytest.raises(VaspSubmitted):
            calc.calculate()

        # Second call - queued
        with pytest.raises(VaspQueued):
            calc.calculate()

        # Third call - running
        with pytest.raises(VaspRunning):
            calc.calculate()

        # Fourth call - complete
        calc.calculate()
        assert 'energy' in calc.results


class TestWorkflowMethods:
    """Test workflow-friendly methods."""

    def test_submit(self, temp_dir, si_atoms, mock_runner):
        """Test submit method."""
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=mock_runner)
        jobid = calc.submit()

        # LocalRunner/MockRunner doesn't return jobid for sync completion
        assert jobid is None or isinstance(jobid, str)

    def test_poll(self, temp_dir, si_atoms, mock_runner):
        """Test poll method."""
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=mock_runner)

        status = calc.poll()
        assert status.state == JobState.NOT_STARTED

    def test_is_complete(self, temp_dir, si_atoms, mock_runner):
        """Test is_complete method."""
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=mock_runner)

        assert not calc.is_complete()

        calc.calculate()
        assert calc.is_complete()

    def test_run_async_complete(self, temp_dir, si_atoms, mock_runner):
        """Test run_async with immediate completion."""
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=mock_runner)

        result = calc.run_async()

        assert isinstance(result, CalculationResult)
        assert result.state == JobState.COMPLETE
        assert result.energy == pytest.approx(-10.84274516)
        assert result.forces is not None

    def test_run_async_submitted(self, temp_dir, si_atoms):
        """Test run_async with submitted state."""
        runner = MockRunner(
            results=MockResults(energy=-10.0),
            state_sequence=[JobState.SUBMITTED]
        )
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=runner)

        result = calc.run_async()

        assert result.state == JobState.SUBMITTED
        assert result.jobid is not None


class TestParameterSet:
    """Test parameter setting."""

    def test_set_single_param(self, vasp_calc):
        """Test setting a single parameter."""
        changed = vasp_calc.set(encut=500)

        assert 'encut' in changed
        assert vasp_calc.parameters['encut'] == 500

    def test_set_clears_results(self, vasp_calc):
        """Test that set() clears results."""
        vasp_calc.calculate()
        assert vasp_calc.results  # Has results

        vasp_calc.set(encut=500)
        assert not vasp_calc.results  # Results cleared

    def test_set_xc_updates_defaults(self, vasp_calc):
        """Test that changing XC updates related parameters."""
        vasp_calc.set(xc='HSE06')

        assert vasp_calc.parameters['lhfcalc'] is True
        assert vasp_calc.parameters['hfscreen'] == 0.2


class TestSerialization:
    """Test pickle serialization for workflow tools."""

    def test_pickle_roundtrip(self, vasp_calc):
        """Test that calculator can be pickled and unpickled."""
        # Calculate to populate results
        vasp_calc.calculate()

        # Pickle and unpickle
        pickled = pickle.dumps(vasp_calc)
        restored = pickle.loads(pickled)

        assert restored.directory == vasp_calc.directory
        assert restored.parameters == vasp_calc.parameters
        assert isinstance(restored.runner, LocalRunner)  # Recreated

    def test_pickle_preserves_results(self, vasp_calc):
        """Test that results are preserved through pickle."""
        vasp_calc.calculate()
        original_energy = vasp_calc.results['energy']

        restored = pickle.loads(pickle.dumps(vasp_calc))

        assert restored.results['energy'] == original_energy


class TestRepr:
    """Test string representation."""

    def test_repr(self, temp_dir, si_atoms):
        """Test __repr__ method."""
        calc = Vasp(label=temp_dir, atoms=si_atoms, xc='HSE06')

        repr_str = repr(calc)

        assert 'Vasp' in repr_str
        assert temp_dir in repr_str
        assert 'HSE06' in repr_str
