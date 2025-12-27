# Calculator

The `Vasp` calculator is the main interface for running VASP calculations.

## Basic Usage

```python
from ase.build import bulk
from vasp import Vasp

atoms = bulk('Si', 'diamond', a=5.43)

calc = Vasp(
    atoms=atoms,
    xc='PBE',
    encut=400,
    kpts=(8, 8, 8),
)

energy = calc.potential_energy
forces = calc.forces
stress = calc.stress
```

## Constructor Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `atoms` | Atoms | ASE Atoms object |
| `label` | str | Calculation directory name |
| `runner` | Runner | Execution backend |
| `xc` | str | Exchange-correlation functional |
| `encut` | float | Plane-wave cutoff (eV) |
| `kpts` | tuple | K-point mesh |
| ... | ... | Any VASP parameter |

## Properties

```python
calc.potential_energy  # Total energy (eV)
calc.forces            # Forces (eV/Å)
calc.stress            # Stress tensor (eV/Å³)
calc.atoms             # Current atoms
calc.parameters        # All parameters
calc.directory         # Working directory
calc.results           # All results dict
```

## Methods

```python
calc.read_results()           # Parse output files
calc.read_doscar()            # Read DOS data
calc.read_procar()            # Read projected DOS
calc.read_eigenval()          # Read eigenvalues
calc.get_work_function()      # Calculate work function
calc.get_band_gap_from_doscar()  # Get band gap
```

## Input File Generation

VASP input files are generated automatically:

- **INCAR**: From keyword parameters
- **POSCAR**: From atoms object
- **KPOINTS**: From `kpts` parameter
- **POTCAR**: From pseudopotential path

## Output Parsing

Results are parsed from VASP output files:

- **OUTCAR**: Energy, forces, stress, magnetic moments
- **vasprun.xml**: All results in structured format
- **DOSCAR**: Density of states
- **EIGENVAL**: Band energies
- **LOCPOT**: Electrostatic potential
