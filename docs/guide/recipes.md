# Recipes

The recipe system provides reusable calculation patterns.

## Core Recipes

### Static Calculation

```python
from vasp.recipes import static_job

result = static_job(atoms, runner=runner, encut=500)

print(result.energy)
print(result.forces)
```

### Relaxation

```python
from vasp.recipes import relax_job

result = relax_job(atoms, runner=runner, relax_cell=True)

relaxed_atoms = result.atoms
```

### Double Relaxation

```python
from vasp.recipes import double_relax_flow

# Coarse → fine relaxation
result = double_relax_flow(atoms, runner=runner)
```

## Slab Recipes

```python
from vasp.recipes.slabs import (
    slab_static_job,
    slab_relax_job,
    bulk_to_slabs_flow,
    make_slabs_from_bulk,
)

# Generate slabs from bulk
slabs = make_slabs_from_bulk(bulk_atoms, miller_indices=[(1,1,1)])

# Calculate all slabs
results = bulk_to_slabs_flow(bulk_atoms, miller_indices=[(1,1,1)])
```

## Phonon Recipes

```python
from vasp.recipes.phonons import phonon_flow

result = phonon_flow(
    atoms,
    runner=runner,
    supercell_matrix=(2, 2, 2),
)

print(result.frequencies)
```

## Custom Recipes

### Create a Job

```python
from vasp.recipes.decorators import job

@job
def my_calculation(atoms, runner=None, **kwargs):
    calc = Vasp(atoms=atoms, runner=runner, **kwargs)
    energy = calc.potential_energy
    return VaspResult(atoms=atoms, energy=energy)
```

### Create a Flow

```python
from vasp.recipes.decorators import flow

@flow
def my_workflow(atoms, runner=None):
    relaxed = relax_job(atoms, runner=runner)
    static = static_job(relaxed.atoms, runner=runner)
    return {'relaxed': relaxed, 'static': static}
```

## Workflow Engines

Set environment variable to enable:

```bash
export VASP_WORKFLOW_ENGINE=prefect
```

Recipes automatically use the engine's decorators.
