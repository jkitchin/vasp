# Migration Guide: v1.x to v2.0

This guide helps you migrate from vasp-ase v1.x (legacy monkeypatch-based) to v2.0 (modern ASE calculator).

## Overview of Changes

Version 2.0 is a **complete rewrite** with significant breaking changes. The new version provides:
- ✅ Modern ASE calculator interface
- ✅ Better type safety with comprehensive type hints
- ✅ Cleaner API without monkeypatching
- ✅ Pluggable execution backends
- ✅ Async workflow support
- ✅ Parameter preset functions
- ✅ Workflow recipes

**Removed features:**
- ❌ MongoDB integration
- ❌ Dynamic keyword validation via `vasp.validate`
- ❌ Configuration via `vasp.vasprc`
- ❌ Legacy getter/setter methods

---

## Installation

### Uninstall v1.x (Optional)

```bash
pip uninstall vasp
```

### Install v2.0

```bash
pip install vasp-ase
```

Or for development:

```bash
git clone https://github.com/jkitchin/vasp.git
cd vasp
pip install -e ".[dev]"
```

---

## Key API Changes

### 1. Calculator Initialization

**Before (v1.x):**
```python
from vasp import Vasp

# Could create calculator without atoms
calc = Vasp('calc-dir')

# Or with atoms
calc = Vasp('calc-dir', xc='PBE', encut=400)

# Set atoms later
calc.set_atoms(atoms)

# Update parameters via set()
calc.set(encut=500, kpts=[8, 8, 8])
```

**After (v2.0):**
```python
from vasp import Vasp

# Atoms required at initialization or set explicitly
calc = Vasp(
    'calc-dir',           # label/directory
    atoms=atoms,          # atoms object
    xc='PBE',
    encut=400
)

# Or set atoms as property
calc = Vasp('calc-dir', xc='PBE')
calc.atoms = atoms

# Parameters passed as keyword arguments
calc = Vasp('calc-dir', atoms=atoms, encut=500, kpts=(8, 8, 8))
```

**Key differences:**
- Directory is now called `label` (first positional argument)
- `atoms` parameter encouraged at initialization
- All parameters as keyword arguments (no `set()` method)
- `kpts` as tuple instead of list

### 2. Running Calculations

**Before (v1.x):**
```python
# Explicit calculate() call
calc.calculate()

# Or get energy
energy = calc.get_potential_energy()
```

**After (v2.0):**
```python
# Access properties directly (triggers calculation if needed)
energy = calc.potential_energy

# Or use ASE's get_ methods (still work)
energy = calc.get_potential_energy()

# No more explicit calculate() - handled automatically
```

### 3. Accessing Results

**Before (v1.x):**
```python
# Various getter methods
energy = calc.get_potential_energy()
forces = calc.get_forces()
atoms = calc.get_atoms()

# Custom getters
fermi = calc.get_fermi_level()
dos = calc.get_dos()
```

**After (v2.0):**
```python
# Standard ASE properties
energy = calc.potential_energy  # or calc.get_potential_energy()
forces = calc.forces            # or calc.get_forces()
atoms = calc.atoms

# Results dict
fermi = calc.results.get('fermi_level')

# Specialized methods from mixins
dos = calc.get_dos(...)
bands = calc.get_band_structure(...)
charges = calc.bader()
```

### 4. Parameter Updates

**Before (v1.x):**
```python
# Update parameters
calc.set(encut=500)
calc.set(kpts=[10, 10, 10], ismear=1)

# Get parameters
params = calc.get_parameters()
```

**After (v2.0):**
```python
# Create new calculator with updated parameters
calc = Vasp(label='new-dir', atoms=atoms, encut=500, kpts=(10, 10, 10))

# Or use set() method (returns changed parameters)
changed = calc.set(encut=500, kpts=(10, 10, 10))

# Access parameters
params = calc.parameters
```

### 5. File I/O

**Before (v1.x):**
```python
# Read existing calculation
calc = Vasp('calc-dir')
calc.read()  # Explicit read

# Write input files
calc.write_input()
```

**After (v2.0):**
```python
# Read existing calculation (automatic)
calc = Vasp('calc-dir')
# Results loaded automatically when accessed

# Write input files (automatic before calculation)
calc.write_input()  # Still available but usually automatic
```

---

## Feature Migration

### Van der Waals Corrections

**Before (v1.x):**
```python
# Manually set vdW parameters
calc = Vasp(
    atoms=atoms,
    xc='PBE',
    ivdw=12,  # D3-BJ
    # ... other params
)
```

**After (v2.0):**
```python
from vasp.parameters import get_vdw_params

# Use preset functions
vdw_params = get_vdw_params('d3bj')
calc = Vasp(atoms=atoms, xc='PBE', **vdw_params)

# Available methods: 'd2', 'd3', 'd3bj', 'ts', 'vdw-df2', 'mbd', etc.
```

### DFT+U

**Before (v1.x):**
```python
# Manual LDAU setup
calc = Vasp(
    atoms=atoms,
    ldau=True,
    ldautype=2,
    ldaul=[2, -1],      # Fe=d, O=none
    ldauu=[4.0, 0.0],
    ldauj=[0.0, 0.0],
)
```

**After (v2.0):**
```python
from vasp.parameters import get_ldau_params, HubbardU

# Use preset with automatic setup
ldau_params = get_ldau_params(
    symbols=['Fe', 'O'],
    u_values={'Fe': HubbardU(u=4.0, j=0.0, l=2)}
)
calc = Vasp(atoms=atoms, **ldau_params)

# Or use simple float (defaults to common values)
ldau_params = get_ldau_params(['Fe', 'O'], {'Fe': 4.0})
```

### Hybrid Functionals

**Before (v1.x):**
```python
# Manual HSE06 setup
calc = Vasp(
    atoms=atoms,
    xc='PBE',
    lhfcalc=True,
    hfscreen=0.2,
    algo='Damped',
    time=0.4,
)
```

**After (v2.0):**
```python
from vasp.parameters import get_hybrid_params

# Use preset
hse_params = get_hybrid_params('hse06')
calc = Vasp(atoms=atoms, **hse_params)

# Available: 'hse06', 'hse03', 'pbe0', 'b3lyp'
```

### Relaxation

**Before (v1.x):**
```python
# Set up relaxation
calc = Vasp(
    atoms=atoms,
    ibrion=2,
    isif=3,
    nsw=100,
    ediffg=-0.02,
)
calc.calculate()
relaxed = calc.get_atoms()
```

**After (v2.0):**
```python
# Option 1: Direct calculator
calc = Vasp(
    atoms=atoms,
    ibrion=2,
    isif=3,
    nsw=100,
    ediffg=-0.02,
)
energy = calc.potential_energy  # Triggers calculation
relaxed = calc.atoms

# Option 2: Use recipe (recommended)
from vasp.recipes import relax_job

result = relax_job(atoms, relax_cell=True)
relaxed = result.atoms
energy = result.energy
```

---

## Removed Features

### MongoDB Integration

**No direct replacement.** MongoDB functionality has been removed.

**Migration strategy:**
1. Export your MongoDB data before upgrading
2. Use the new vector database for structure similarity
3. Or implement custom storage using the calculator results

```python
# v2.0 - Custom storage example
import json

calc = Vasp('calc', atoms=atoms, xc='PBE')
energy = calc.potential_energy

# Save results manually
results = {
    'directory': calc.directory,
    'atoms': atoms.todict(),
    'parameters': calc.parameters,
    'energy': energy,
    'forces': calc.forces.tolist(),
}

with open('results.json', 'w') as f:
    json.dump(results, f)
```

**New alternative - Vector Database:**
```python
from vasp.database import VectorDatabase, CrystalEmbedder

db = VectorDatabase()
embedder = CrystalEmbedder()

# Store calculation with embedding
db.add_structure(atoms, calc=calc, embedder=embedder)

# Find similar structures
matches = db.find_similar_structures(query_atoms, k=5, embedder=embedder)
```

### Configuration via vasprc

**No direct replacement.** Configuration is now done via:

1. **Environment variables:**
   ```bash
   export VASP_PP_PATH=/path/to/potpaw_PBE
   export ASE_VASP_COMMAND="mpirun -np 4 vasp_std"
   ```

2. **Runner objects:**
   ```python
   from vasp.runners import LocalRunner

   runner = LocalRunner(
       command='mpirun -np 8 vasp_std',
       pp_path='/custom/potcar/path'
   )
   calc = Vasp(atoms=atoms, runner=runner, ...)
   ```

### Legacy Validation

**No direct replacement.** Parameter validation has been simplified.

**Before (v1.x):**
```python
from vasp.validate import keywords, keyword_alist

# Get all valid keywords
all_keywords = keywords()

# Get documentation
docs = keyword_alist()
```

**After (v2.0):**
- Type hints provide static validation
- Runtime validation happens in calculator
- Use `/vasp-help` command for parameter documentation
- Check VASP wiki: https://www.vasp.at/wiki/

---

## New Features to Adopt

### 1. Runners for Different Environments

```python
from vasp.runners import LocalRunner, SlurmRunner, MockRunner

# Local execution
local_runner = LocalRunner(command='mpirun -np 4 vasp_std')

# SLURM cluster
slurm_runner = SlurmRunner(
    partition='compute',
    nodes=2,
    ntasks_per_node=24,
    time='2:00:00'
)

# Testing without VASP
from vasp.runners import MockResults
mock_runner = MockRunner(results=MockResults(energy=-10.5))

# Use with calculator
calc = Vasp(atoms=atoms, runner=slurm_runner, ...)
```

### 2. Workflow Recipes

```python
from vasp.recipes import static_job, relax_job, double_relax_flow

# Single-point energy
result = static_job(atoms)

# Geometry optimization
result = relax_job(atoms, relax_cell=True)

# Production relaxation (coarse then fine)
result = double_relax_flow(atoms)

# Access results
print(result.energy)
print(result.atoms)
```

### 3. Async Workflow Pattern

```python
from vasp.exceptions import VaspSubmitted, VaspQueued, VaspRunning

calc = Vasp('calc', atoms=atoms, runner=slurm_runner, ...)

try:
    energy = calc.potential_energy
except VaspSubmitted as e:
    print(f"Job {e.jobid} submitted - continue with other work")
except VaspQueued:
    print("Job in queue - check back later")
except VaspRunning:
    print("Job is running - check back later")
```

### 4. Claude Code Integration

```bash
# Install global skills
pip install vasp-ase
vasp-claude install

# Use commands
/vasp-help encut           # Parameter help
/watch-job calc-dir        # Monitor job
/fix-job calc-dir          # Diagnose errors
```

---

## Common Migration Patterns

### Pattern 1: Simple Energy Calculation

**Before:**
```python
from vasp import Vasp
from ase.build import bulk

atoms = bulk('Si')
calc = Vasp('si-calc', xc='PBE', encut=400, kpts=[8,8,8])
calc.set_atoms(atoms)
calc.calculate()
energy = calc.get_potential_energy()
```

**After:**
```python
from vasp import Vasp
from ase.build import bulk

atoms = bulk('Si')
calc = Vasp('si-calc', atoms=atoms, xc='PBE', encut=400, kpts=(8,8,8))
energy = calc.potential_energy
```

### Pattern 2: Geometry Optimization

**Before:**
```python
calc = Vasp('relax', xc='PBE', ibrion=2, isif=3, nsw=100, ediffg=-0.02)
calc.set_atoms(atoms)
calc.calculate()
relaxed = calc.get_atoms()
```

**After:**
```python
# Option 1: Direct
calc = Vasp('relax', atoms=atoms, xc='PBE',
            ibrion=2, isif=3, nsw=100, ediffg=-0.02)
energy = calc.potential_energy
relaxed = calc.atoms

# Option 2: Recipe (recommended)
from vasp.recipes import relax_job
result = relax_job(atoms, relax_cell=True, directory='relax')
relaxed = result.atoms
```

### Pattern 3: DOS Calculation

**Before:**
```python
# SCF
calc1 = Vasp('scf', xc='PBE', ...)
calc1.set_atoms(atoms)
calc1.calculate()

# DOS
calc2 = Vasp('dos', xc='PBE', ismear=-5, nedos=2000, icharg=11)
calc2.set_atoms(atoms)
calc2.calculate()
dos = calc2.get_dos()
```

**After:**
```python
# SCF
calc_scf = Vasp('scf', atoms=atoms, xc='PBE', ...)
energy = calc_scf.potential_energy

# DOS
calc_dos = Vasp('dos', atoms=atoms, xc='PBE',
                ismear=-5, nedos=2000, icharg=11)
dos = calc_dos.get_dos()
```

### Pattern 4: Reading Existing Calculations

**Before:**
```python
calc = Vasp('existing-calc')
calc.read()
energy = calc.get_potential_energy()
```

**After:**
```python
calc = Vasp('existing-calc')
# Results loaded automatically
energy = calc.potential_energy
```

---

## Troubleshooting

### Issue: "ModuleNotFoundError: No module named 'vasp.validate'"

**Solution:** The `validate` module was removed. Remove imports and use:
- Parameter presets: `from vasp.parameters import get_vdw_params, ...`
- VASP documentation: https://www.vasp.at/wiki/
- Claude help: `/vasp-help <parameter>`

### Issue: "AttributeError: 'Vasp' object has no attribute 'set'"

**Solution:** Pass parameters directly to constructor or use `calc.set(...)` which is still available:
```python
# Old way still works
calc.set(encut=500)

# But creating new calculator is preferred
calc = Vasp(label='new', atoms=atoms, encut=500, ...)
```

### Issue: "ModuleNotFoundError: No module named 'vasp.mongo'"

**Solution:** MongoDB integration was removed. Use:
- Vector database: `from vasp.database import VectorDatabase`
- Custom JSON/pickle storage
- External database of your choice

### Issue: Parameters not recognized

**Solution:** Check parameter names (may have changed):
- `kpts=[8,8,8]` → `kpts=(8,8,8)` (tuple instead of list)
- All parameter names are lowercase
- Check type hints for expected types

### Issue: "Calculator must have atoms"

**Solution:** Always provide atoms at initialization or set explicitly:
```python
# At initialization
calc = Vasp('dir', atoms=atoms, ...)

# Or set property
calc = Vasp('dir', ...)
calc.atoms = atoms
```

---

## Testing Your Migration

Create a test script to verify basic functionality:

```python
"""Test v2.0 migration."""
from ase.build import bulk
from vasp import Vasp
from vasp.parameters import get_vdw_params
from vasp.runners import MockRunner, MockResults

# Test 1: Basic calculation
atoms = bulk('Si')
mock = MockRunner(results=MockResults(energy=-10.5))
calc = Vasp('test', atoms=atoms, runner=mock, xc='PBE', encut=400)
energy = calc.potential_energy
assert energy == -10.5
print("✓ Basic calculation works")

# Test 2: Parameter presets
vdw = get_vdw_params('d3bj')
calc = Vasp('test2', atoms=atoms, runner=mock, **vdw)
assert 'ivdw' in calc.parameters
print("✓ Parameter presets work")

# Test 3: Recipes
from vasp.recipes import static_job
result = static_job(atoms, runner=mock, directory='test3')
assert result.energy == -10.5
print("✓ Recipes work")

print("\n✅ Migration successful!")
```

---

## Getting Help

- **Documentation**: https://github.com/jkitchin/vasp
- **Examples**: Check `examples/` directory (21 progressive tutorials)
- **Issues**: https://github.com/jkitchin/vasp/issues
- **Claude Code**: Use `/vasp-help <topic>` after running `vasp-claude install`

---

## Summary Checklist

- [ ] Update `from vasp import Vasp` imports
- [ ] Add `atoms` parameter to Vasp() initialization
- [ ] Replace `calc.calculate()` with property access (`calc.potential_energy`)
- [ ] Change `kpts=[...]` to `kpts=(...)`
- [ ] Replace `calc.set_atoms()` with `atoms=` parameter
- [ ] Remove MongoDB dependencies
- [ ] Remove `vasp.validate` imports
- [ ] Remove `vasp.vasprc` usage
- [ ] Consider using parameter presets (`get_vdw_params`, etc.)
- [ ] Consider using workflow recipes for common tasks
- [ ] Update tests to use MockRunner
- [ ] Test migrated code

Welcome to vasp-ase 2.0! 🎉
