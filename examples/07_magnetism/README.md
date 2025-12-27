# 07 - Magnetism

Calculate magnetic properties of iron and nickel oxide.

## Key Concepts

- **Spin Polarization**: ISPIN=2 for magnetic calculations
- **Magnetic Moment**: Net spin moment in Bohr magnetons (μB)
- **Ferromagnetism**: Parallel spin alignment (Fe)
- **Antiferromagnetism**: Alternating spin alignment (NiO)
- **DFT+U**: Correction for correlated d/f electrons

## VASP Parameters

### Ferromagnetic Fe

```python
calc = Vasp(
    atoms=fe,
    xc='PBE',
    encut=400,
    kpts=(12, 12, 12),
    ismear=1,           # Metal smearing
    sigma=0.1,
    ispin=2,            # Spin-polarized
    magmom=[3.0],       # Initial moment
    lorbit=11,          # Orbital decomposition
)
```

### NiO with DFT+U

```python
from vasp.parameters import get_ldau_params, HubbardU

ldau = get_ldau_params(['Ni', 'O'], {'Ni': HubbardU(u=6.45)})

calc = Vasp(
    atoms=nio,
    xc='PBE',
    encut=500,
    kpts=(8, 8, 8),
    ispin=2,
    magmom=[2.0, 0.0],
    **ldau,
)
```

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `ispin=2` | Spin-polarized | Enable spin polarization |
| `magmom` | [M1, M2, ...] | Initial magnetic moments |
| `lorbit=11` | Full decomposition | Atom/orbital moments |
| `ldauu` | [U1, U2, ...] | Hubbard U values |
| `ldaul` | [l1, l2, ...] | Angular momentum for U |

## Physical Background

### Ferromagnetic Iron
- BCC crystal structure
- 3d⁶4s² electron configuration
- Magnetic moment: ~2.2 μB/atom
- Curie temperature: 1043 K

### Antiferromagnetic NiO
- Rock salt structure
- Type-II AFM ordering along [111]
- Ni moment: ~1.9 μB
- Néel temperature: 523 K

### DFT+U Correction
Standard DFT underestimates d-electron localization. DFT+U adds:
- Hubbard U: On-site Coulomb repulsion
- Hund's J: Exchange parameter

Typical U values:
- Ni (oxides): 6-7 eV
- Fe (oxides): 4-5 eV
- Co (oxides): 3-4 eV

## Output

```
Fe (BCC):
  Magnetic moment: 2.18 μB (exp: 2.22 μB)
  Magnetic stabilization: 0.92 eV

NiO:
  DFT magnetic moment: 1.45 μB
  DFT+U magnetic moment: 1.72 μB
```

## Tips

1. **Initial Moments**: Set reasonable starting values
2. **Convergence**: Magnetic calculations may need more iterations
3. **U Values**: Use literature values or compute with linear response
4. **AFM Ordering**: Requires appropriate supercell and magnetic structure

## Applications

- Magnetic metals and alloys
- Transition metal oxides
- Magnetic anisotropy
- Exchange coupling constants

## Next Steps

Try `08_surfaces/` for slab and surface calculations.
