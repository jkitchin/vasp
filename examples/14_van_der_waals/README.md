# 14 - Van der Waals Corrections

Apply dispersion corrections for systems with London dispersion forces using graphite.

## Key Concepts

- **Dispersion Forces**: Weak attraction from induced dipoles
- **D3-BJ**: Grimme's D3 with Becke-Johnson damping
- **vdW-DF**: Nonlocal density functional for dispersion
- **Layered Materials**: Where vdW is essential

## VASP Parameters

### D3-BJ Correction

```python
from vasp.parameters import get_vdw_params

d3bj = get_vdw_params('d3bj')

calc = Vasp(
    atoms=graphite,
    xc='PBE',
    encut=500,
    kpts=(12, 12, 6),
    ismear=0,
    sigma=0.05,
    **d3bj,
)
```

### Available Methods

| Method | IVDW | Description |
|--------|------|-------------|
| D2 | 10 | Grimme D2, fixed C6 |
| D3 | 11 | D3 with zero-damping |
| D3-BJ | 12 | D3 with BJ damping (recommended) |
| TS | 20 | Tkatchenko-Scheffler |
| TS-SCS | 202 | TS with self-consistent screening |
| dDsC | 4 | density-dependent dispersion |

## Physical Background

### Why DFT Fails

Standard DFT (LDA/GGA):
- Based on local density approximation
- Cannot describe long-range correlation (London dispersion)
- Gives ~0 binding for graphite layers (wrong!)

### Dispersion Corrections

**D3-BJ formula:**
```
E_disp = -Σ_ij [C6_ij/(R6_ij + f(R0)) + C8_ij/(R8_ij + f(R0))]
```

Features:
- Geometry-dependent C6 coefficients
- Becke-Johnson damping at short range
- No double-counting at covalent distances

## Output

```
Graphite interlayer binding:

Method      Interlayer (Å)    Binding (meV/atom)
PBE         Unbound           ~0
PBE-D3(BJ)  3.35              ~50
vdW-DF2     3.40              ~48
Experiment  3.35              ~52
```

## Tips

1. **Default Choice**: D3-BJ for most applications
2. **Metals**: TS-SCS often better for metallic systems
3. **Cost**: D2/D3 nearly free, TS/vdW-DF more expensive
4. **Consistency**: Use same method for all calculations in study

## When to Use

**Essential for:**
- Layered materials (graphite, MoS2, etc.)
- Molecular crystals
- Adsorption on surfaces
- Biomolecular systems
- Noble gas interactions

**Not needed for:**
- Covalent/ionic bonding only
- Bulk metals (but helps for surfaces)
- When dispersion is small compared to other effects

## Common Systems

| System | Correction Needed? |
|--------|-------------------|
| Graphite | Yes, essential |
| MoS2 | Yes |
| CO on Pt | Yes, for binding energy |
| Bulk Si | No |
| Water on surface | Yes |

## Applications

- Interlayer binding in 2D materials
- Molecular adsorption energies
- Crystal structure prediction
- Protein-ligand binding

## Next Steps

Try `15_workflows/` for automated calculation workflows.
