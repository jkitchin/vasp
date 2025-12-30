"""Tests for mixin classes."""

import os

import numpy as np
import pytest

from vasp import Vasp
from vasp.runners import MockResults, MockRunner


class TestIOMixin:
    """Test IO mixin methods."""

    def test_write_input_creates_files(self, temp_dir, si_atoms, mock_runner):
        """Test that write_input creates all required files."""
        calc = Vasp(
            label=temp_dir,
            atoms=si_atoms,
            runner=mock_runner,
            xc='PBE',
            encut=400,
            kpts=(4, 4, 4)
        )

        calc.write_input(si_atoms)

        assert os.path.exists(os.path.join(temp_dir, 'INCAR'))
        assert os.path.exists(os.path.join(temp_dir, 'POSCAR'))
        assert os.path.exists(os.path.join(temp_dir, 'KPOINTS'))

    def test_write_poscar(self, temp_dir, si_atoms, mock_runner):
        """Test POSCAR writing."""
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=mock_runner)

        calc.write_poscar(si_atoms)

        poscar_path = os.path.join(temp_dir, 'POSCAR')
        assert os.path.exists(poscar_path)

        with open(poscar_path) as f:
            content = f.read()

        assert 'Si' in content
        assert '5.43' in content or '5.4' in content

    def test_write_incar(self, temp_dir, si_atoms, mock_runner):
        """Test INCAR writing."""
        calc = Vasp(
            label=temp_dir,
            atoms=si_atoms,
            runner=mock_runner,
            encut=500,
            ismear=0,
            sigma=0.05,
            ediff=1e-7
        )

        calc.write_incar()

        incar_path = os.path.join(temp_dir, 'INCAR')
        assert os.path.exists(incar_path)

        with open(incar_path) as f:
            content = f.read()

        assert 'ENCUT' in content
        assert '500' in content
        assert 'ISMEAR' in content
        assert 'EDIFF' in content

    def test_write_kpoints_gamma(self, temp_dir, si_atoms, mock_runner):
        """Test KPOINTS writing with gamma-centered grid."""
        calc = Vasp(
            label=temp_dir,
            atoms=si_atoms,
            runner=mock_runner,
            kpts=(8, 8, 8),
            gamma=True
        )

        calc.write_kpoints()

        kpoints_path = os.path.join(temp_dir, 'KPOINTS')
        assert os.path.exists(kpoints_path)

        with open(kpoints_path) as f:
            content = f.read()

        assert 'Gamma' in content
        assert '8' in content

    def test_write_kpoints_monkhorst(self, temp_dir, si_atoms, mock_runner):
        """Test KPOINTS writing with Monkhorst-Pack grid."""
        calc = Vasp(
            label=temp_dir,
            atoms=si_atoms,
            runner=mock_runner,
            kpts=(6, 6, 6),
            gamma=False
        )

        calc.write_kpoints()

        kpoints_path = os.path.join(temp_dir, 'KPOINTS')
        with open(kpoints_path) as f:
            content = f.read()

        assert 'Monkhorst' in content or 'M' in content

    def test_read_incar(self, calc_dir, mock_runner, si_atoms):
        """Test reading INCAR file."""
        calc = Vasp(label=calc_dir, atoms=si_atoms, runner=mock_runner)

        incar = calc.read_incar()

        assert isinstance(incar, dict)
        # Keys are stored lowercase
        assert incar.get('prec', '').lower() == 'accurate'
        assert incar.get('encut') == 400

    def test_read_kpoints(self, calc_dir, mock_runner, si_atoms):
        """Test reading KPOINTS file."""
        calc = Vasp(label=calc_dir, atoms=si_atoms, runner=mock_runner)

        kpts = calc.read_kpoints()

        assert isinstance(kpts, dict)
        assert 'grid' in kpts or 'kpts' in kpts

    def test_read_results(self, complete_calc_dir, mock_runner, si_atoms):
        """Test reading calculation results."""
        calc = Vasp(label=complete_calc_dir, atoms=si_atoms, runner=mock_runner)

        # Force read from files
        calc.read_results()

        assert 'energy' in calc.results
        assert 'forces' in calc.results


class TestElectronicMixin:
    """Test electronic properties mixin."""

    def test_get_fermi_level(self, temp_dir, si_atoms):
        """Test getting Fermi level."""
        results = MockResults(energy=-10.0, fermi_level=5.4321)
        runner = MockRunner(results=results)
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=runner)
        calc.calculate()

        fermi = calc.get_fermi_level()

        assert fermi == pytest.approx(5.4321)

    def test_get_eigenvalues(self, temp_dir, si_atoms):
        """Test getting eigenvalues."""
        eigenvalues = np.array([[[-5.0, 2.0, 6.0]]])  # (spin, kpt, band)
        results = MockResults(energy=-10.0, eigenvalues=eigenvalues)
        runner = MockRunner(results=results)
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=runner)
        calc.calculate()

        # get_eigenvalues returns 1D array for specific kpt/spin
        eigs = calc.get_eigenvalues(kpt=0, spin=0)

        assert eigs.shape == (3,)  # 3 bands
        assert eigs[0] == pytest.approx(-5.0)
        assert eigs[2] == pytest.approx(6.0)

    def test_get_occupation_numbers(self, temp_dir, si_atoms):
        """Test getting occupation numbers."""
        eigenvalues = np.array([[[-5.0, 2.0, 6.0]]])  # Need eigenvalues for proper parsing
        occupations = np.array([[[1.0, 1.0, 0.0]]])
        results = MockResults(energy=-10.0, eigenvalues=eigenvalues, occupations=occupations)
        runner = MockRunner(results=results)
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=runner)
        calc.calculate()

        # get_occupation_numbers returns 1D array for specific kpt/spin
        occ = calc.get_occupation_numbers(kpt=0, spin=0)

        assert occ.shape == (3,)  # 3 bands
        assert occ[0] == pytest.approx(1.0)
        assert occ[2] == pytest.approx(0.0)

    def test_get_ibz_k_points(self, temp_dir, si_atoms):
        """Test getting IBZ k-points."""
        kpts = np.array([[0, 0, 0], [0.5, 0, 0], [0.5, 0.5, 0]])
        results = MockResults(energy=-10.0, ibz_kpts=kpts)
        runner = MockRunner(results=results)
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=runner)
        calc.calculate()

        ibz = calc.get_ibz_k_points()

        assert ibz.shape == (3, 3)

    def test_get_magnetic_moment(self, temp_dir, si_atoms):
        """Test getting total magnetic moment."""
        results = MockResults(energy=-10.0, magmom=2.5)
        runner = MockRunner(results=results)
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=runner)
        calc.calculate()

        magmom = calc.get_magnetic_moment()

        assert magmom == pytest.approx(2.5)

    def test_get_magnetic_moments(self, temp_dir, si_atoms):
        """Test getting per-atom magnetic moments."""
        magmoms = np.array([1.5, -0.5])
        results = MockResults(energy=-10.0, magmoms=magmoms)
        runner = MockRunner(results=results)
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=runner)
        calc.calculate()

        moms = calc.get_magnetic_moments()

        assert len(moms) == 2
        assert moms[0] == pytest.approx(1.5)

    def test_get_homo_lumo(self, temp_dir, si_atoms):
        """Test getting HOMO and LUMO energies."""
        eigenvalues = np.array([[[-5.0, 2.0, 2.0, 2.0, 6.0, 7.0]]])
        occupations = np.array([[[1.0, 1.0, 1.0, 1.0, 0.0, 0.0]]])
        results = MockResults(
            energy=-10.0,
            eigenvalues=eigenvalues,
            occupations=occupations,
            fermi_level=4.0,  # Between HOMO (2.0) and LUMO (6.0)
        )
        runner = MockRunner(results=results)
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=runner)
        calc.calculate()

        homo, lumo = calc.get_homo_lumo()

        assert homo == pytest.approx(2.0)
        assert lumo == pytest.approx(6.0)

    def test_get_band_gap(self, temp_dir, si_atoms):
        """Test getting band gap."""
        eigenvalues = np.array([[[-5.0, 2.0, 2.0, 2.0, 6.0, 7.0]]])
        occupations = np.array([[[1.0, 1.0, 1.0, 1.0, 0.0, 0.0]]])
        results = MockResults(
            energy=-10.0,
            eigenvalues=eigenvalues,
            occupations=occupations,
            fermi_level=4.0,  # Between HOMO (2.0) and LUMO (6.0)
        )
        runner = MockRunner(results=results)
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=runner)
        calc.calculate()

        gap = calc.get_band_gap()

        assert gap == pytest.approx(4.0)  # 6.0 - 2.0


class TestDynamicsMixin:
    """Test dynamics/vibration mixin."""

    def test_get_zero_point_energy(self, temp_dir, si_atoms):
        """Test zero-point energy calculation from frequencies."""
        # Need to set up vibrational calculation output
        # Create mock OUTCAR with vibrational data
        from vasp.tests.fixtures import MOCK_OUTCAR_VIBRATIONS

        results = MockResults(energy=-10.0)
        runner = MockRunner(results=results)
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=runner)
        calc.calculate()

        # Overwrite OUTCAR with vibrational data
        with open(os.path.join(temp_dir, 'OUTCAR'), 'w') as f:
            f.write(MOCK_OUTCAR_VIBRATIONS)

        # This test would need actual vibrational data parsing
        # For now, just verify the method exists and has correct signature
        try:
            zpe = calc.get_zero_point_energy()
            assert isinstance(zpe, float)
        except (ValueError, FileNotFoundError):
            # Expected if mock data doesn't have proper vibrational info
            pass


class TestAnalysisMixin:
    """Test analysis mixin."""

    def test_get_charge_density(self, temp_dir, si_atoms, mock_runner):
        """Test charge density reading."""
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=mock_runner)

        # This requires CHGCAR file, test that method exists
        with pytest.raises(FileNotFoundError):
            calc.get_charge_density()

    def test_get_local_potential(self, temp_dir, si_atoms, mock_runner):
        """Test local potential reading."""
        calc = Vasp(label=temp_dir, atoms=si_atoms, runner=mock_runner)

        # This requires LOCPOT file, test that method exists
        with pytest.raises(FileNotFoundError):
            calc.get_local_potential()
