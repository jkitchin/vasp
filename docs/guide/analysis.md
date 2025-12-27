# Analysis

Methods for analyzing VASP calculation results.

## Reading Output Files

### DOSCAR

```python
dos_data = calc.read_doscar()

energy = dos_data['energy']
total_dos = dos_data['total_dos']
fermi = dos_data['fermi']
projected = dos_data.get('projected')  # If LORBIT=11
```

### EIGENVAL

```python
eigen_data = calc.read_eigenval()

eigenvalues = eigen_data['eigenvalues']
kpoints = eigen_data['kpoints']
weights = eigen_data['weights']
```

### PROCAR

```python
procar_data = calc.read_procar()

# Orbital projections for band structure coloring
```

### LOCPOT

```python
# Work function from electrostatic potential
work_function = calc.get_work_function()
```

## Band Gap Analysis

```python
# From DOSCAR
band_gap = calc.get_band_gap_from_doscar(tol=0.1)

# From EIGENVAL
eigenvalues = calc.read_eigenval()['eigenvalues']
fermi = calc.results['fermi_level']

occupied = eigenvalues[eigenvalues < fermi]
unoccupied = eigenvalues[eigenvalues > fermi]
vbm = occupied.max()
cbm = unoccupied.min()
gap = cbm - vbm
```

## Magnetic Moments

```python
# Total magnetic moment
total_mag = calc.results.get('magnetic_moment')

# Per-atom moments (with LORBIT=11)
atom_mags = calc.results.get('magnetic_moments')
```

## Surface Analysis

```python
from vasp.recipes.slabs import calculate_surface_energy

# Calculate surface energy
e_surf = calculate_surface_energy(slab_result, bulk_energy_per_atom)
```

## Equation of State

```python
from ase.eos import EquationOfState

volumes = [...]
energies = [...]

eos = EquationOfState(volumes, energies, eos='birchmurnaghan')
v0, e0, B = eos.fit()

print(f"Equilibrium volume: {v0:.2f} Å³")
print(f"Bulk modulus: {B * 160.2:.1f} GPa")
```
