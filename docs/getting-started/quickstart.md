# Quick Start

## Basic Calculation

```python
from ase.build import bulk
from vasp import Vasp

# Create structure
atoms = bulk('Si', 'diamond', a=5.43)

# Create calculator
calc = Vasp(
    atoms=atoms,
    xc='PBE',
    encut=400,
    kpts=(8, 8, 8),
)

# Get energy
energy = calc.potential_energy
print(f"Energy: {energy:.3f} eV")
```

## Structure Relaxation

```python
calc = Vasp(
    atoms=atoms,
    xc='PBE',
    encut=400,
    kpts=(8, 8, 8),
    ibrion=2,      # Conjugate gradient
    isif=3,        # Relax cell and positions
    nsw=100,       # Max ionic steps
    ediffg=-0.01,  # Force convergence
)

energy = calc.potential_energy
relaxed_atoms = calc.atoms
```

## Using Recipes

```python
from vasp.recipes import relax_job, static_job

# Relax structure
result = relax_job(atoms, relax_cell=True)

# Static calculation on relaxed structure
static = static_job(result.atoms, lorbit=11)

print(f"Relaxed energy: {result.energy:.3f} eV")
```

## Using Parameter Presets

```python
from vasp.parameters import get_vdw_params, get_ldau_params, HubbardU

# Add D3-BJ dispersion
vdw = get_vdw_params('d3bj')
calc = Vasp(atoms=atoms, **vdw)

# DFT+U for transition metal oxides
ldau = get_ldau_params(['Fe', 'O'], {'Fe': HubbardU(u=4.0)})
calc = Vasp(atoms=atoms, **ldau)
```

## HPC Execution

```python
from vasp import Vasp
from vasp.runners import SlurmRunner

runner = SlurmRunner(
    partition='compute',
    nodes=2,
    time='4:00:00',
)

calc = Vasp(atoms=atoms, runner=runner, ...)

try:
    energy = calc.potential_energy
except VaspSubmitted as e:
    print(f"Job submitted: {e.jobid}")
```

## Next Steps

See the [Tutorials](../tutorials/beginner.md) for detailed examples.
