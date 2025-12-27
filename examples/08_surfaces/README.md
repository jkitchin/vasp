# 08 - Surfaces

Calculate surface properties of Cu(111), including surface energy and work function.

## Key Concepts

- **Slab Model**: Periodic slab with vacuum separation
- **Surface Energy**: Energy cost to create surface
- **Work Function**: Minimum energy to remove electron
- **Dipole Correction**: Corrects for asymmetric slabs
- **Layer Convergence**: Surface properties vs. slab thickness

## VASP Parameters

### Slab Calculation

```python
from ase.build import fcc111
from ase.constraints import FixAtoms

slab = fcc111('Cu', size=(2, 2, 4), vacuum=12.0)

# Fix bottom 2 layers
constraint = FixAtoms(mask=positions_z <= z_threshold)
slab.set_constraint(constraint)

calc = Vasp(
    atoms=slab,
    xc='PBE',
    encut=400,
    kpts=(6, 6, 1),     # 1 k-point in z
    ismear=1,
    sigma=0.1,
    ldipol=True,        # Dipole correction
    idipol=3,           # Correct in z
    isif=2,             # Relax positions only
    ibrion=2,
    nsw=50,
    ediffg=-0.02,
)
```

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `vacuum` | 12-15 Å | Separation between images |
| `ldipol=True` | Dipole correction | For asymmetric slabs |
| `idipol=3` | z-direction | Correct along surface normal |
| `isif=2` | Fix cell | Only relax atomic positions |
| `lvtot=True` | Write LOCPOT | For work function |

## Physical Background

### Surface Energy

Surface energy γ is the energy cost to create surface area:

```
γ = (E_slab - N × E_bulk) / (2 × A)
```

where:
- E_slab: Total slab energy
- N: Number of atoms in slab
- E_bulk: Energy per bulk atom
- A: Surface area (factor 2 for both sides)

### Work Function

Work function Φ is:

```
Φ = V_vacuum - E_Fermi
```

Measured from the planar-averaged electrostatic potential in the vacuum region.

## Output

```
Cu(111) Surface Properties:
  Surface energy: 1.25 J/m² (exp: ~1.83 J/m²)
  Work function: ~4.9 eV (exp: 4.94 eV)
```

Note: GGA tends to underestimate surface energies.

## Convergence Tests

Always check convergence with:
1. **Slab thickness**: 4-6 layers typically sufficient
2. **Vacuum thickness**: 10-15 Å to avoid image interaction
3. **K-points**: Dense in-plane grid
4. **Fixed layers**: Fix bottom 2-3 layers

## Tips

1. **Symmetric vs. Asymmetric**: Use symmetric slabs when possible
2. **Relaxation**: Only top layers should relax
3. **Dipole Correction**: Essential for polar surfaces
4. **K-grid**: Use Γ-centered with 1 in z-direction

## Applications

- Surface stability comparison
- Surface reconstruction
- Work function engineering
- Catalysis studies

## Next Steps

Try `09_adsorption/` for molecules on surfaces.
