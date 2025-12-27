---
name: vasp-calculation
description: Create a new VASP calculation using the vasp-ase interface. Use when the user wants to set up a DFT calculation, compute energies, or run VASP.
---

# VASP Calculation Skill

When helping users create VASP calculations:

## 1. Understand the Task

Ask clarifying questions if needed:
- What property to calculate (energy, relaxation, DOS, bands, etc.)?
- What material/structure?
- What level of theory (PBE, HSE06, with vdW, etc.)?

## 2. Create the Structure

Use ASE to create or load the structure:

```python
from ase.build import bulk, molecule, fcc111
from ase.io import read

# For crystals
atoms = bulk('Si', 'diamond', a=5.43)

# For molecules
atoms = molecule('H2O')
atoms.center(vacuum=10.0)

# For surfaces
slab = fcc111('Cu', size=(2, 2, 4), vacuum=10.0)

# From file
atoms = read('structure.cif')
```

## 3. Set Up the Calculator

Use appropriate parameters:

```python
from vasp import Vasp
from vasp.parameters import get_vdw_params, get_ldau_params, get_hybrid_params

# Basic calculation
calc = Vasp(
    atoms=atoms,
    xc='PBE',
    encut=400,
    kpts=(8, 8, 8),
    ismear=1,
    sigma=0.1,
)

# With vdW correction
vdw = get_vdw_params('d3bj')
calc = Vasp(atoms=atoms, xc='PBE', encut=400, **vdw)

# With DFT+U
ldau = get_ldau_params(['Fe', 'O'], {'Fe': HubbardU(u=4.0)})
calc = Vasp(atoms=atoms, xc='PBE', **ldau)
```

## 4. Run and Get Results

```python
# Energy
energy = calc.potential_energy

# Forces
forces = calc.forces

# Stress
stress = calc.stress

# For relaxation
calc = Vasp(
    atoms=atoms,
    ibrion=2,
    isif=3,
    nsw=100,
    ediffg=-0.02,
)
energy = calc.potential_energy
relaxed_atoms = calc.atoms
```

## 5. Use Recipes for Common Workflows

```python
from vasp.recipes import static_job, relax_job, double_relax_flow

# Quick static calculation
result = static_job(atoms)

# Relaxation
result = relax_job(atoms, relax_cell=True)

# Production-quality relaxation
result = double_relax_flow(atoms)
```

## Key Parameters Reference

| Task | Key Parameters |
|------|---------------|
| Static | ENCUT, KPTS, ISMEAR, SIGMA |
| Relax | IBRION=2, ISIF=3, NSW, EDIFFG |
| DOS | ISMEAR=-5, NEDOS, LORBIT=11, ICHARG=11 |
| Bands | ICHARG=11, dense k-path |
| Magnetic | ISPIN=2, MAGMOM |
| Hybrid | LHFCALC=True, HFSCREEN, ALGO=Damped |
