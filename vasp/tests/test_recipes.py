"""Tests for recipe module."""

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk, fcc111

from vasp.recipes.core import VaspResult, double_relax_flow, relax_job, static_job
from vasp.recipes.decorators import WORKFLOW_ENGINE, flow, job, subflow
from vasp.recipes.phonons import PhononResult, _check_phonopy, get_phonopy_supercells
from vasp.recipes.slabs import (
    SlabResult,
    bulk_to_slabs_flow,
    make_slabs_from_bulk,
    slab_relax_job,
    slab_static_job,
)
from vasp.runners import MockResults, MockRunner


class TestDecorators:
    """Test recipe decorators."""

    def test_job_decorator_no_engine(self):
        """Test @job decorator without workflow engine."""
        @job
        def my_job(x: int) -> int:
            return x * 2

        result = my_job(5)
        assert result == 10

    def test_flow_decorator_no_engine(self):
        """Test @flow decorator without workflow engine."""
        @flow
        def my_flow(x: int) -> int:
            return x + 10

        result = my_flow(5)
        assert result == 15

    def test_subflow_decorator_no_engine(self):
        """Test @subflow decorator without workflow engine."""
        @subflow
        def my_subflow(items: list) -> list:
            return [x * 2 for x in items]

        result = my_subflow([1, 2, 3])
        assert result == [2, 4, 6]

    def test_job_with_kwargs(self):
        """Test @job with keyword arguments."""
        @job(name='custom_job')
        def custom_job(x: int) -> int:
            return x * 3

        result = custom_job(4)
        assert result == 12

    def test_workflow_engine_detection(self):
        """Test that workflow engine is detected from environment."""
        # WORKFLOW_ENGINE is set at import time from env var
        # Valid values are None or one of the supported engines
        # Allow any value since it comes from environment
        assert WORKFLOW_ENGINE is None or isinstance(WORKFLOW_ENGINE, str)


class TestVaspResult:
    """Test VaspResult dataclass."""

    def test_vasp_result_creation(self):
        """Test creating VaspResult."""
        atoms = bulk('Si')
        result = VaspResult(
            atoms=atoms,
            energy=-10.5,
            forces=np.zeros((2, 3)),
            stress=np.zeros(6),
        )
        assert result.energy == -10.5
        assert len(result.atoms) == 2

    def test_vasp_result_defaults(self):
        """Test VaspResult default values."""
        result = VaspResult()
        assert result.atoms is None
        assert result.energy is None
        assert result.forces is None
        assert result.converged is True


class TestCoreRecipes:
    """Test core VASP recipes."""

    @pytest.fixture
    def si_atoms(self):
        """Create silicon unit cell."""
        return bulk('Si', 'diamond', a=5.43)

    @pytest.fixture
    def mock_runner(self):
        """Create mock runner for testing."""
        results = MockResults(
            energy=-10.5,
            forces=np.zeros((2, 3)),
            stress=np.zeros(6),
        )
        return MockRunner(results=results)

    def test_static_job(self, si_atoms, mock_runner):
        """Test static calculation job."""
        result = static_job(si_atoms, runner=mock_runner)

        assert isinstance(result, VaspResult)
        assert result.energy == pytest.approx(-10.5)
        assert result.atoms is not None
        assert len(result.atoms) == len(si_atoms)

    def test_static_job_with_params(self, si_atoms, mock_runner):
        """Test static job with custom parameters."""
        result = static_job(
            si_atoms,
            runner=mock_runner,
            encut=500,
            kpts=(8, 8, 8),
        )

        assert result.energy is not None
        assert result.parameters.get('encut') == 500

    def test_relax_job(self, si_atoms, mock_runner):
        """Test relaxation job."""
        result = relax_job(si_atoms, runner=mock_runner)

        assert isinstance(result, VaspResult)
        assert result.energy is not None

    def test_relax_job_forces(self, si_atoms, mock_runner):
        """Test relaxation job returns forces."""
        result = relax_job(si_atoms, runner=mock_runner)
        assert result.forces is not None

    def test_double_relax_flow(self, si_atoms, mock_runner):
        """Test double relaxation flow."""
        result = double_relax_flow(si_atoms, runner=mock_runner)

        assert isinstance(result, VaspResult)
        assert result.energy is not None


class TestSlabResult:
    """Test SlabResult dataclass."""

    def test_slab_result_creation(self):
        """Test creating SlabResult."""
        atoms = fcc111('Cu', size=(2, 2, 4), vacuum=10.0)
        result = SlabResult(
            atoms=atoms,
            energy=-50.0,
            miller_indices=(1, 1, 1),
            surface_energy=0.05,
        )
        assert result.energy == -50.0
        assert result.miller_indices == (1, 1, 1)
        assert result.surface_energy == 0.05


class TestSlabRecipes:
    """Test slab VASP recipes."""

    @pytest.fixture
    def cu_slab(self):
        """Create copper slab."""
        return fcc111('Cu', size=(2, 2, 4), vacuum=10.0)

    @pytest.fixture
    def cu_bulk(self):
        """Create copper bulk."""
        return bulk('Cu', 'fcc', a=3.6)

    @pytest.fixture
    def mock_runner(self, tmp_path):
        """Create mock runner for testing."""
        results = MockResults(
            energy=-50.0,
            forces=np.random.randn(16, 3) * 0.01,
            stress=np.zeros(6),
        )
        return MockRunner(results=results)

    def test_slab_static_job(self, cu_slab, mock_runner):
        """Test slab static calculation."""
        result = slab_static_job(cu_slab, runner=mock_runner)

        assert isinstance(result, VaspResult)
        assert result.energy is not None

    def test_slab_relax_job(self, cu_slab, mock_runner):
        """Test slab relaxation job."""
        result = slab_relax_job(cu_slab, runner=mock_runner)

        assert isinstance(result, VaspResult)
        assert result.energy is not None

    def test_make_slabs_from_bulk(self, cu_bulk):
        """Test slab generation from bulk."""
        slabs = make_slabs_from_bulk(
            cu_bulk,
            miller_indices=[(1, 0, 0), (1, 1, 0), (1, 1, 1)],
            min_slab_size=8.0,
            min_vacuum_size=10.0,
        )

        assert len(slabs) >= 1
        for slab in slabs:
            assert isinstance(slab, Atoms)
            assert slab.pbc[2]  # Should have vacuum in z (pbc=True)

    def test_bulk_to_slabs_flow(self, cu_bulk, mock_runner):
        """Test bulk to slabs workflow."""
        # This flow generates slabs and runs calculations
        results = bulk_to_slabs_flow(
            cu_bulk,
            runner=mock_runner,
            miller_indices=[(1, 1, 1)],
        )

        # Returns dict with 'bulk' and 'slabs' keys
        assert isinstance(results, dict)
        assert 'bulk' in results
        assert 'slabs' in results


class TestPhononResult:
    """Test PhononResult dataclass."""

    def test_phonon_result_creation(self):
        """Test creating PhononResult."""
        result = PhononResult(
            supercell_matrix=(2, 2, 2),
            frequencies=[100.0, 200.0, 300.0],
            directory='/tmp/phonon',
        )
        assert result.supercell_matrix == (2, 2, 2)
        assert len(result.frequencies) == 3


class TestPhononRecipes:
    """Test phonon recipes."""

    @pytest.fixture
    def si_atoms(self):
        """Create silicon unit cell."""
        return bulk('Si', 'diamond', a=5.43)

    def test_check_phonopy(self):
        """Test phonopy availability check."""
        result = _check_phonopy()
        # Result is boolean
        assert isinstance(result, bool)

    @pytest.mark.skipif(not _check_phonopy(), reason="phonopy not installed")
    def test_get_phonopy_supercells(self, si_atoms):
        """Test generating phonopy supercells."""
        phonon, displaced_cells = get_phonopy_supercells(
            si_atoms,
            supercell_matrix=(2, 2, 2),
            displacement=0.01,
        )

        assert len(displaced_cells) > 0
        # Each displaced cell should be larger than unit cell
        for cell in displaced_cells:
            assert len(cell) > len(si_atoms)

    @pytest.mark.skipif(not _check_phonopy(), reason="phonopy not installed")
    def test_phonopy_supercell_sizes(self, si_atoms):
        """Test supercell matrix affects size correctly."""
        _, cells_2x2x2 = get_phonopy_supercells(si_atoms, supercell_matrix=(2, 2, 2))
        _, cells_3x3x3 = get_phonopy_supercells(si_atoms, supercell_matrix=(3, 3, 3))

        # 3x3x3 cells should be larger
        assert len(cells_3x3x3[0]) > len(cells_2x2x2[0])


class TestRecipeIntegration:
    """Integration tests for recipes."""

    @pytest.fixture
    def si_atoms(self):
        """Create silicon unit cell."""
        return bulk('Si', 'diamond', a=5.43)

    @pytest.fixture
    def mock_runner(self, tmp_path):
        """Create mock runner."""
        results = MockResults(
            energy=-10.5,
            forces=np.zeros((2, 3)),
            stress=np.zeros(6),
        )
        return MockRunner(results=results)

    def test_recipe_chain(self, si_atoms, mock_runner):
        """Test chaining multiple recipes."""
        # First relax
        relax_result = relax_job(si_atoms, runner=mock_runner)

        # Then static on relaxed structure
        static_result = static_job(relax_result.atoms, runner=mock_runner)

        assert relax_result.energy is not None
        assert static_result.energy is not None

    def test_recipe_with_vdw(self, si_atoms, mock_runner):
        """Test recipe with VdW corrections."""
        from vasp.parameters import get_vdw_params

        vdw_params = get_vdw_params('d3bj')
        result = static_job(si_atoms, runner=mock_runner, **vdw_params)

        assert result.energy is not None
        # VdW parameters should be in calc parameters
        assert result.parameters.get('ivdw') == 12

    def test_recipe_with_ldau(self, si_atoms, mock_runner):
        """Test recipe with DFT+U (conceptual - Si doesn't need U)."""
        from vasp.parameters import HubbardU, get_ldau_params

        # This is just for testing the parameter passing
        # Si doesn't actually need U corrections
        # get_ldau_params(symbols, u_values) - symbols first, then u_values
        ldau_params = get_ldau_params(['Si'], {'Si': HubbardU(u=0.0)})
        result = static_job(si_atoms, runner=mock_runner, **ldau_params)

        assert result.energy is not None
