# 18 - 3D Visualization

Visualize volumetric data from VASP: charge density, electrostatic potential, and electron localization function (ELF).

## Key Concepts

- **Charge Density (CHGCAR)**: Electron density distribution ρ(r)
- **Local Potential (LOCPOT)**: Electrostatic potential V(r)
- **ELF (ELFCAR)**: Electron localization function (0 to 1)
- **Isosurfaces**: 3D contour surfaces at constant value

## VASP Output Files

| File | Content | Settings Required |
|------|---------|-------------------|
| `CHGCAR` | Total charge density | Default (LCHARG=.TRUE.) |
| `LOCPOT` | Local potential | LVTOT=.TRUE. or LVHAR=.TRUE. |
| `ELFCAR` | Electron localization | LELF=.TRUE. |
| `PARCHG` | Partial charge | LPARD=.TRUE. + band selection |

## VASP Parameters

```python
calc = Vasp(
    atoms=atoms,
    xc='PBE',
    encut=400,
    kpts=(8, 8, 8),
    lcharg=True,    # Write CHGCAR (default)
    lvtot=True,     # Write LOCPOT (total potential)
    lelf=True,      # Write ELFCAR (electron localization)
    prec='Accurate',
)
```

### Key Settings

| Parameter | Value | Description |
|-----------|-------|-------------|
| `lcharg` | True | Write charge density |
| `lvtot` | True | Write local potential (ionic + Hartree + XC) |
| `lvhar` | True | Write Hartree potential only |
| `lelf` | True | Write electron localization function |
| `ngx/ngy/ngz` | Auto | FFT grid (controls resolution) |

## Physical Background

### Charge Density

Total electron density from occupied Kohn-Sham orbitals:

```
ρ(r) = Σₙ fₙ |ψₙ(r)|²
```

Used to analyze bonding, charge transfer, Bader analysis.

### Electrostatic Potential

From Poisson equation:

```
∇²V(r) = -4πρ(r)
```

Used for work function, band alignment, interface dipoles.

### Electron Localization Function (ELF)

Measures electron pairing:

```
ELF = 1 / (1 + (D/D_h)²)
```

- ELF = 1: Perfect localization (lone pairs, covalent bonds)
- ELF = 0.5: Electron gas-like
- ELF → 0: Very delocalized

## Visualization Tools

### Python (matplotlib + ASE)

```python
from ase.io.cube import read_cube
import matplotlib.pyplot as plt

# Read VASP volumetric data
data, atoms = read_cube('CHGCAR')

# 2D slice
plt.contourf(data[:, :, nz//2])
```

### Python (py3Dmol for 3D)

```python
import py3Dmol

view = py3Dmol.view()
view.addVolumetricData(data, "cube", {'isoval': 0.1})
view.show()
```

### External Tools

- **VESTA**: Crystal structure viewer with isosurface support
- **ParaView**: Scientific visualization, large datasets
- **VMD**: Molecular dynamics visualization
- **Avogadro**: Molecule editor with visualization

## Output Examples

### Charge Density Isosurface

```
Si charge density at ρ = 0.05 e/Å³:
  - Electron density concentrated at bond centers
  - Typical covalent bonding pattern
  - Valence electrons between Si atoms
```

### Work Function from LOCPOT

```
Vacuum level: 4.32 eV (from LOCPOT plateau)
Fermi level: -0.85 eV
Work function: 4.32 - (-0.85) = 5.17 eV
Experimental Cu(111): 4.94 eV
```

### ELF Analysis

```
ELF values in H2O:
  - O lone pairs: ELF ≈ 0.95
  - O-H bonds: ELF ≈ 0.98
  - Between atoms: ELF ≈ 0.0
```

## Tips

1. **Grid Resolution**: Increase NGX/NGY/NGZ for smoother isosurfaces
2. **Isosurface Value**: Experiment to find informative values
3. **Large Files**: CHGCAR can be >100 MB; use compression or sampling
4. **Spin**: CHGCAR contains both spin channels if ISPIN=2
5. **Units**: VASP uses electrons/Å³ for density, eV for potential

## Common Visualizations

| Property | File | Typical Isovalues |
|----------|------|-------------------|
| Valence density | CHGCAR | 0.02-0.1 e/Å³ |
| Bond charge | CHG | 0.05-0.2 e/Å³ |
| ELF | ELFCAR | 0.7-0.9 |
| Work function | LOCPOT | Slice through surface |

## File Format Notes

VASP uses a special format for volumetric data:
- Header: lattice vectors, atom positions
- Data: values on FFT grid (fortran ordering)

ASE can read/write these files with `ase.io.vasp`.

## Applications

- Chemical bonding visualization
- Charge transfer analysis
- Work function calculation
- Interface dipole analysis
- Bader charge analysis (with bader code)
- Reaction mechanism visualization

## Next Steps

Congratulations! You've completed the tutorial series. Explore:
- Combining techniques for complex materials problems
- Automating workflows with `vasp.recipes`
- Using the Claude Code integration for assistance
