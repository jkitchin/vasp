# 15 - Workflows

Automated calculation workflows using the recipe system.

## Key Concepts

- **Recipes**: Reusable calculation patterns
- **@job**: Single VASP calculation
- **@flow**: Multi-step workflow
- **Workflow Engines**: Integration with Prefect, Dask, etc.

## Recipe System

### Basic Usage

```python
from vasp.recipes import static_job, relax_job, double_relax_flow

# Single-point calculation
result = static_job(atoms, runner=runner)

# Geometry optimization
result = relax_job(atoms, runner=runner, relax_cell=True)

# Coarse → fine relaxation
result = double_relax_flow(atoms, runner=runner)
```

### Available Recipes

| Recipe | Type | Description |
|--------|------|-------------|
| `static_job` | @job | Single-point energy |
| `relax_job` | @job | Geometry optimization |
| `double_relax_flow` | @flow | Coarse → fine relax |
| `slab_static_job` | @job | Surface static calc |
| `slab_relax_job` | @job | Surface relaxation |
| `bulk_to_slabs_flow` | @flow | Full slab workflow |
| `phonon_flow` | @flow | Phonon calculation |

## Creating Custom Recipes

### Custom Job

```python
from vasp.recipes.decorators import job

@job
def my_custom_job(atoms, runner=None, **kwargs):
    """Custom calculation pattern."""
    calc = Vasp(atoms=atoms, runner=runner, **kwargs)
    energy = calc.potential_energy
    return VaspResult(
        atoms=atoms.copy(),
        energy=energy,
        parameters=calc.parameters.copy(),
    )
```

### Custom Flow

```python
from vasp.recipes.decorators import flow

@flow
def my_analysis_flow(atoms, runner=None):
    """Multi-step analysis."""
    # Step 1
    relaxed = relax_job(atoms, runner=runner)

    # Step 2
    dos = static_job(relaxed.atoms, runner=runner, lorbit=11)

    return {'relaxed': relaxed, 'dos': dos}
```

## Workflow Engine Integration

Set environment variable to enable:

```bash
export VASP_WORKFLOW_ENGINE=prefect  # or dask, parsl, covalent
```

### Supported Engines

| Engine | Use Case |
|--------|----------|
| Prefect | Cloud workflows, monitoring |
| Dask | Distributed computing |
| Parsl | HPC systems |
| Covalent | Cloud/hybrid execution |
| Jobflow | Materials Project integration |

### Prefect Example

```python
# With VASP_WORKFLOW_ENGINE=prefect

from vasp.recipes import relax_job, static_job

# These are now Prefect tasks/flows!
# Features: retries, caching, logging, web UI
```

## Result Dataclasses

### VaspResult

```python
@dataclass
class VaspResult:
    atoms: Atoms | None = None
    energy: float | None = None
    forces: np.ndarray | None = None
    stress: np.ndarray | None = None
    fermi_level: float | None = None
    band_gap: float | None = None
    magnetic_moment: float | None = None
    parameters: dict = field(default_factory=dict)
    directory: str | None = None
    converged: bool = True
```

### SlabResult

```python
@dataclass
class SlabResult(VaspResult):
    miller_indices: tuple[int, int, int] | None = None
    surface_energy: float | None = None
    work_function: float | None = None
    layers: int = 0
```

## Testing with MockRunner

```python
from vasp.runners import MockRunner, MockResults

# Create mock results
mock = MockResults(energy=-10.5, forces=np.zeros((2, 3)))
runner = MockRunner(results=mock)

# Test without VASP
result = static_job(atoms, runner=runner)
```

## Best Practices

1. **Use Recipes**: More reproducible than ad-hoc scripts
2. **Set Runners**: Configure for your HPC environment
3. **MockRunner**: Test workflows without VASP
4. **Workflow Engine**: Enable for production HPC use
5. **Custom Recipes**: Create for repeated patterns

## Production Example

```python
from vasp.recipes import double_relax_flow, static_job
from vasp.runners import SlurmRunner

runner = SlurmRunner(partition='compute', nodes=2)

def screen_material(atoms):
    # Relax
    relaxed = double_relax_flow(atoms, runner=runner)

    # Properties
    props = static_job(
        relaxed.atoms,
        runner=runner,
        lorbit=11,
        nedos=2000,
    )

    return {
        'structure': relaxed.atoms,
        'energy': relaxed.energy,
        'band_gap': props.band_gap,
    }
```

## Applications

- High-throughput screening
- Automated convergence testing
- Complex multi-step calculations
- Reproducible research workflows

## Conclusion

You've completed all 15 examples! Key takeaways:

1. Start simple (01-05)
2. Add complexity as needed (06-10)
3. Use advanced features appropriately (11-15)
4. Leverage the recipe system for reproducibility
