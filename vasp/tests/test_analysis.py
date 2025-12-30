"""Tests for analysis mixin methods."""

import os

import pytest

from vasp import Vasp
from vasp.runners import MockResults, MockRunner

# Mock file contents
MOCK_DOSCAR = """     2     2     1     1
  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00
  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00
  0.0000000000000000E+00  0.0000000000000000E+00
    1.0000000     0.0000000
  10.00000  -10.00000   100   5.5000   1.0000
 -10.0000   0.0500    0.00500
  -9.8000   0.1000    0.01500
  -9.6000   0.2000    0.03500
  -9.4000   0.3000    0.06500
  -9.2000   0.5000    0.11500
  -9.0000   0.8000    0.19500
  -8.8000   1.2000    0.31500
  -8.6000   1.8000    0.49500
  -8.4000   2.5000    0.74500
  -8.2000   3.2000    1.06500
  -8.0000   4.0000    1.46500
  -7.8000   4.5000    1.91500
  -7.6000   5.0000    2.41500
  -7.4000   5.2000    2.93500
  -7.2000   5.0000    3.43500
  -7.0000   4.5000    3.88500
  -6.8000   4.0000    4.28500
  -6.6000   3.5000    4.63500
  -6.4000   3.0000    4.93500
  -6.2000   2.5000    5.18500
  -6.0000   2.0000    5.38500
  -5.8000   1.5000    5.53500
  -5.6000   1.0000    5.63500
  -5.4000   0.6000    5.69500
  -5.2000   0.3000    5.72500
  -5.0000   0.1000    5.73500
  -4.8000   0.0500    5.74000
  -4.6000   0.0200    5.74200
  -4.4000   0.0100    5.74300
  -4.2000   0.0050    5.74350
  -4.0000   0.0020    5.74370
  -3.8000   0.0010    5.74380
  -3.6000   0.0005    5.74385
  -3.4000   0.0002    5.74387
  -3.2000   0.0001    5.74388
  -3.0000   0.0000    5.74388
  -2.8000   0.0000    5.74388
  -2.6000   0.0000    5.74388
  -2.4000   0.0000    5.74388
  -2.2000   0.0000    5.74388
  -2.0000   0.0000    5.74388
  -1.8000   0.0000    5.74388
  -1.6000   0.0000    5.74388
  -1.4000   0.0000    5.74388
  -1.2000   0.0000    5.74388
  -1.0000   0.0000    5.74388
  -0.8000   0.0000    5.74388
  -0.6000   0.0000    5.74388
  -0.4000   0.0000    5.74388
  -0.2000   0.0000    5.74388
   0.0000   0.0000    5.74388
   0.2000   0.0000    5.74388
   0.4000   0.0000    5.74388
   0.6000   0.0000    5.74388
   0.8000   0.0000    5.74388
   1.0000   0.0000    5.74388
   1.2000   0.0001    5.74389
   1.4000   0.0002    5.74391
   1.6000   0.0005    5.74396
   1.8000   0.0010    5.74406
   2.0000   0.0020    5.74426
   2.2000   0.0050    5.74476
   2.4000   0.0100    5.74576
   2.6000   0.0200    5.74776
   2.8000   0.0500    5.75276
   3.0000   0.1000    5.76276
   3.2000   0.2000    5.78276
   3.4000   0.3000    5.81276
   3.6000   0.5000    5.86276
   3.8000   0.8000    5.94276
   4.0000   1.2000    6.06276
   4.2000   1.8000    6.24276
   4.4000   2.5000    6.49276
   4.6000   3.2000    6.81276
   4.8000   4.0000    7.21276
   5.0000   4.5000    7.66276
   5.2000   5.0000    8.16276
   5.4000   5.2000    8.68276
   5.6000   5.0000    9.18276
   5.8000   4.5000    9.63276
   6.0000   4.0000   10.03276
   6.2000   3.5000   10.38276
   6.4000   3.0000   10.68276
   6.6000   2.5000   10.93276
   6.8000   2.0000   11.13276
   7.0000   1.5000   11.28276
   7.2000   1.0000   11.38276
   7.4000   0.6000   11.44276
   7.6000   0.3000   11.47276
   7.8000   0.1000   11.48276
   8.0000   0.0500   11.48776
   8.2000   0.0200   11.48976
   8.4000   0.0100   11.49076
   8.6000   0.0050   11.49126
   8.8000   0.0020   11.49146
   9.0000   0.0010   11.49156
   9.2000   0.0005   11.49161
   9.4000   0.0002   11.49163
   9.6000   0.0001   11.49164
   9.8000   0.0000   11.49164
  10.0000   0.0000   11.49164
"""

MOCK_PROCAR = """PROCAR lm decomposed + phase
# of k-points:    3         # of bands:    4         # of ions:     2

 k-point     1 :    0.00000000 0.00000000 0.00000000     weight = 0.12500000

band     1 # energy   -6.52430000 # occ.  2.00000000

ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot
  1  0.050  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.050
  2  0.050  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.050
tot  0.100  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.100

band     2 # energy   -2.10230000 # occ.  2.00000000

ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot
  1  0.030  0.010  0.010  0.010  0.000  0.000  0.000  0.000  0.000  0.060
  2  0.030  0.010  0.010  0.010  0.000  0.000  0.000  0.000  0.000  0.060
tot  0.060  0.020  0.020  0.020  0.000  0.000  0.000  0.000  0.000  0.120

band     3 # energy    3.45670000 # occ.  0.00000000

ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot
  1  0.010  0.020  0.020  0.020  0.000  0.000  0.000  0.000  0.000  0.070
  2  0.010  0.020  0.020  0.020  0.000  0.000  0.000  0.000  0.000  0.070
tot  0.020  0.040  0.040  0.040  0.000  0.000  0.000  0.000  0.000  0.140

band     4 # energy    5.67890000 # occ.  0.00000000

ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot
  1  0.005  0.030  0.030  0.030  0.000  0.000  0.000  0.000  0.000  0.095
  2  0.005  0.030  0.030  0.030  0.000  0.000  0.000  0.000  0.000  0.095
tot  0.010  0.060  0.060  0.060  0.000  0.000  0.000  0.000  0.000  0.190

 k-point     2 :    0.50000000 0.00000000 0.00000000     weight = 0.37500000

band     1 # energy   -5.12340000 # occ.  2.00000000

ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot
  1  0.045  0.005  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.050
  2  0.045  0.005  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.050
tot  0.090  0.010  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.100

band     2 # energy   -1.23450000 # occ.  2.00000000

ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot
  1  0.025  0.015  0.010  0.010  0.000  0.000  0.000  0.000  0.000  0.060
  2  0.025  0.015  0.010  0.010  0.000  0.000  0.000  0.000  0.000  0.060
tot  0.050  0.030  0.020  0.020  0.000  0.000  0.000  0.000  0.000  0.120

band     3 # energy    4.56780000 # occ.  0.00000000

ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot
  1  0.008  0.025  0.020  0.020  0.000  0.000  0.000  0.000  0.000  0.073
  2  0.008  0.025  0.020  0.020  0.000  0.000  0.000  0.000  0.000  0.073
tot  0.016  0.050  0.040  0.040  0.000  0.000  0.000  0.000  0.000  0.146

band     4 # energy    6.78900000 # occ.  0.00000000

ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot
  1  0.003  0.035  0.030  0.030  0.000  0.000  0.000  0.000  0.000  0.098
  2  0.003  0.035  0.030  0.030  0.000  0.000  0.000  0.000  0.000  0.098
tot  0.006  0.070  0.060  0.060  0.000  0.000  0.000  0.000  0.000  0.196

 k-point     3 :    0.50000000 0.50000000 0.00000000     weight = 0.50000000

band     1 # energy   -4.56780000 # occ.  2.00000000

ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot
  1  0.040  0.005  0.005  0.000  0.000  0.000  0.000  0.000  0.000  0.050
  2  0.040  0.005  0.005  0.000  0.000  0.000  0.000  0.000  0.000  0.050
tot  0.080  0.010  0.010  0.000  0.000  0.000  0.000  0.000  0.000  0.100

band     2 # energy   -0.12340000 # occ.  2.00000000

ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot
  1  0.020  0.015  0.015  0.010  0.000  0.000  0.000  0.000  0.000  0.060
  2  0.020  0.015  0.015  0.010  0.000  0.000  0.000  0.000  0.000  0.060
tot  0.040  0.030  0.030  0.020  0.000  0.000  0.000  0.000  0.000  0.120

band     3 # energy    5.67890000 # occ.  0.00000000

ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot
  1  0.006  0.028  0.022  0.020  0.000  0.000  0.000  0.000  0.000  0.076
  2  0.006  0.028  0.022  0.020  0.000  0.000  0.000  0.000  0.000  0.076
tot  0.012  0.056  0.044  0.040  0.000  0.000  0.000  0.000  0.000  0.152

band     4 # energy    7.89010000 # occ.  0.00000000

ion      s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot
  1  0.002  0.038  0.032  0.030  0.000  0.000  0.000  0.000  0.000  0.102
  2  0.002  0.038  0.032  0.030  0.000  0.000  0.000  0.000  0.000  0.102
tot  0.004  0.076  0.064  0.060  0.000  0.000  0.000  0.000  0.000  0.204
"""


def generate_mock_locpot():
    """Generate a mock LOCPOT file content."""
    header = """Slab calculation
   1.00000000000000
     5.0000000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000    5.0000000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000   20.0000000000000000
   Cu
     4
Direct
  0.00  0.00  0.25
  0.50  0.50  0.25
  0.00  0.50  0.35
  0.50  0.00  0.35

   4   4  20
"""
    # Generate potential values (vacuum at edges, bulk in middle)
    values = []
    for z in range(20):
        for _y in range(4):
            for _x in range(4):
                # Slab centered around z=10
                # Vacuum potential is higher
                if z < 4 or z > 16:
                    val = 5.0  # Vacuum level
                else:
                    val = 0.5  # Bulk region
                values.append(val)

    # Format values
    lines = []
    for i in range(0, len(values), 5):
        line_vals = values[i:i + 5]
        lines.append("  ".join(f"{v:12.8f}" for v in line_vals))

    return header + "\n".join(lines)


@pytest.fixture
def si_atoms():
    """Create silicon unit cell."""
    from ase.build import bulk
    return bulk('Si', 'diamond', a=5.43)


class TestDOSCARReading:
    """Test DOSCAR parsing."""

    def test_read_doscar(self, si_atoms, tmp_path):
        """Test reading DOSCAR file."""
        results = MockResults(energy=-10.5, fermi_level=5.5)
        runner = MockRunner(results=results)
        calc = Vasp(label=str(tmp_path), atoms=si_atoms, runner=runner)

        # Write mock DOSCAR
        doscar_path = os.path.join(str(tmp_path), 'DOSCAR')
        with open(doscar_path, 'w') as f:
            f.write(MOCK_DOSCAR)

        dos_data = calc.read_doscar()

        assert 'energy' in dos_data
        assert 'total_dos' in dos_data
        assert 'fermi' in dos_data
        assert dos_data['fermi'] == pytest.approx(5.5)
        assert len(dos_data['energy']) == 100

    def test_read_doscar_energy_range(self, si_atoms, tmp_path):
        """Test DOSCAR energy range."""
        runner = MockRunner()
        calc = Vasp(label=str(tmp_path), atoms=si_atoms, runner=runner)

        doscar_path = os.path.join(str(tmp_path), 'DOSCAR')
        with open(doscar_path, 'w') as f:
            f.write(MOCK_DOSCAR)

        dos_data = calc.read_doscar()

        # Check energy range (100 points from -10 to ~10)
        assert dos_data['energy'][0] == pytest.approx(-10.0)
        assert dos_data['energy'][-1] == pytest.approx(10.0, abs=0.5)  # Last point may not be exactly 10.0

    def test_read_doscar_not_found(self, si_atoms, tmp_path):
        """Test error when DOSCAR not found."""
        runner = MockRunner()
        calc = Vasp(label=str(tmp_path), atoms=si_atoms, runner=runner)

        with pytest.raises(FileNotFoundError, match="DOSCAR not found"):
            calc.read_doscar()


class TestPROCARReading:
    """Test PROCAR parsing."""

    def test_read_procar(self, si_atoms, tmp_path):
        """Test reading PROCAR file."""
        runner = MockRunner()
        calc = Vasp(label=str(tmp_path), atoms=si_atoms, runner=runner)

        # Write mock PROCAR
        procar_path = os.path.join(str(tmp_path), 'PROCAR')
        with open(procar_path, 'w') as f:
            f.write(MOCK_PROCAR)

        procar_data = calc.read_procar()

        assert procar_data['nkpts'] == 3
        assert procar_data['nbands'] == 4
        assert procar_data['nions'] == 2
        assert len(procar_data['kpoints']) == 3
        assert len(procar_data['weights']) == 3

    def test_read_procar_kpoints(self, si_atoms, tmp_path):
        """Test PROCAR k-point parsing."""
        runner = MockRunner()
        calc = Vasp(label=str(tmp_path), atoms=si_atoms, runner=runner)

        procar_path = os.path.join(str(tmp_path), 'PROCAR')
        with open(procar_path, 'w') as f:
            f.write(MOCK_PROCAR)

        procar_data = calc.read_procar()

        # Check first k-point is gamma
        assert procar_data['kpoints'][0, 0] == pytest.approx(0.0)
        assert procar_data['kpoints'][0, 1] == pytest.approx(0.0)
        assert procar_data['kpoints'][0, 2] == pytest.approx(0.0)

    def test_read_procar_energies(self, si_atoms, tmp_path):
        """Test PROCAR band energy parsing."""
        runner = MockRunner()
        calc = Vasp(label=str(tmp_path), atoms=si_atoms, runner=runner)

        procar_path = os.path.join(str(tmp_path), 'PROCAR')
        with open(procar_path, 'w') as f:
            f.write(MOCK_PROCAR)

        procar_data = calc.read_procar()

        # Check energies shape
        assert procar_data['energies'].shape == (3, 4)
        # First band at gamma
        assert procar_data['energies'][0, 0] == pytest.approx(-6.5243)

    def test_read_procar_not_found(self, si_atoms, tmp_path):
        """Test error when PROCAR not found."""
        runner = MockRunner()
        calc = Vasp(label=str(tmp_path), atoms=si_atoms, runner=runner)

        with pytest.raises(FileNotFoundError, match="PROCAR not found"):
            calc.read_procar()


class TestBandGapFromDOS:
    """Test band gap calculation from DOSCAR."""

    def test_get_band_gap_semiconductor(self, si_atoms, tmp_path):
        """Test band gap for semiconductor."""
        runner = MockRunner()
        calc = Vasp(label=str(tmp_path), atoms=si_atoms, runner=runner)

        doscar_path = os.path.join(str(tmp_path), 'DOSCAR')
        with open(doscar_path, 'w') as f:
            f.write(MOCK_DOSCAR)

        gap = calc.get_band_gap_from_doscar()

        # Silicon has a gap (our mock has a gap around Fermi)
        # The DOS is zero from about -3 eV to +1 eV relative to Fermi
        assert gap > 0
        assert gap < 10.0  # Reasonable range

    def test_get_band_gap_metal(self, si_atoms, tmp_path):
        """Test band gap for metal (should be small/zero)."""
        runner = MockRunner()
        calc = Vasp(label=str(tmp_path), atoms=si_atoms, runner=runner)

        # Create metallic DOSCAR (non-zero DOS at Fermi level)
        # Use fine energy grid to ensure Fermi level is covered
        metal_doscar = """     1     1     1     1
  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00
  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00
  0.0000000000000000E+00  0.0000000000000000E+00
    1.0000000     0.0000000
  10.00000  -10.00000   201   0.0000   1.0000
"""
        # Add DOS data that has states at Fermi level
        # Fine grid with 201 points from -10 to 10 (step=0.1)
        for i in range(201):
            e = -10.0 + i * 0.1
            dos = 1.0  # Non-zero DOS everywhere
            int_dos = i * 0.1
            metal_doscar += f"  {e:.4f}   {dos:.4f}   {int_dos:.4f}\n"

        doscar_path = os.path.join(str(tmp_path), 'DOSCAR')
        with open(doscar_path, 'w') as f:
            f.write(metal_doscar)

        gap = calc.get_band_gap_from_doscar()

        # Metal should have zero or very small gap (within grid resolution)
        assert gap < 0.5  # Small gap due to discretization


class TestWorkFunction:
    """Test work function calculation."""

    def test_get_work_function(self, tmp_path):
        """Test work function calculation from LOCPOT."""
        from ase.build import fcc111

        slab = fcc111('Cu', size=(1, 1, 4), vacuum=10.0)
        results = MockResults(energy=-50.0, fermi_level=4.5)
        runner = MockRunner(results=results)
        calc = Vasp(label=str(tmp_path), atoms=slab, runner=runner)
        calc.calculate()

        # Write mock LOCPOT
        locpot_path = os.path.join(str(tmp_path), 'LOCPOT')
        with open(locpot_path, 'w') as f:
            f.write(generate_mock_locpot())

        # Calculate work function
        wf = calc.get_work_function()

        # Just check that it returns a reasonable value
        # The actual value depends on the mock LOCPOT details
        assert isinstance(wf, float)
        # Work function should be a reasonable value (typically 2-6 eV for metals)
        # Our mock may give different values depending on averaging
        assert abs(wf) < 10.0  # Sanity check


class TestVolumetricFiles:
    """Test reading volumetric data files."""

    def test_get_local_potential(self, tmp_path):
        """Test reading LOCPOT."""
        from ase.build import fcc111

        slab = fcc111('Cu', size=(1, 1, 4), vacuum=10.0)
        runner = MockRunner()
        calc = Vasp(label=str(tmp_path), atoms=slab, runner=runner)

        locpot_path = os.path.join(str(tmp_path), 'LOCPOT')
        with open(locpot_path, 'w') as f:
            f.write(generate_mock_locpot())

        grid, potential = calc.get_local_potential()

        assert grid.shape[0] == 3  # x, y, z grids
        assert potential.shape == (4, 4, 20)

    def test_get_elf_not_found(self, si_atoms, tmp_path):
        """Test error when ELFCAR not found."""
        runner = MockRunner()
        calc = Vasp(label=str(tmp_path), atoms=si_atoms, runner=runner)

        with pytest.raises(FileNotFoundError, match="ELFCAR not found"):
            calc.get_elf()

    def test_get_charge_density_not_found(self, si_atoms, tmp_path):
        """Test error when CHGCAR not found."""
        runner = MockRunner()
        calc = Vasp(label=str(tmp_path), atoms=si_atoms, runner=runner)

        with pytest.raises(FileNotFoundError, match="CHGCAR not found"):
            calc.get_charge_density()
