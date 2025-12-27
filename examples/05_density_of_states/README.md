# 05 - Density of States

Calculate and plot the electronic density of states (DOS) for silicon.

## Key Concepts

- **Density of States**: Distribution of electronic states vs. energy
- **Projected DOS (PDOS)**: Contribution from specific atoms/orbitals
- **Band Gap**: Energy gap between valence and conduction bands
- **Two-Step Calculation**: SCF → non-SCF DOS

## VASP Parameters

### Step 1: Self-Consistent Calculation

```python
calc_scf = Vasp(
    atoms=atoms,
    xc='PBE',
    encut=400,
    kpts=(8, 8, 8),
    ismear=1,
    sigma=0.1,
    lwave=True,   # Save wavefunctions
    lcharg=True,  # Save charge density
)
```

### Step 2: Non-Self-Consistent DOS

```python
calc_dos = Vasp(
    atoms=atoms,
    xc='PBE',
    encut=400,
    kpts=(12, 12, 12),  # Denser k-grid
    ismear=-5,          # Tetrahedron method
    icharg=11,          # Read charge from file
    lorbit=11,          # Atom-projected DOS
    nedos=2001,         # Fine energy grid
)
```

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `ismear=-5` | Tetrahedron | Best for DOS (requires ≥4 k-points) |
| `icharg=11` | Read charge | Non-self-consistent calculation |
| `lorbit=11` | Full PDOS | Atom and orbital decomposition |
| `nedos` | 2001 | Number of energy points |

## Output Files

- **DOSCAR**: Total and projected DOS data
- **EIGENVAL**: Band energies at each k-point

## Physical Background

The DOS g(E) counts states at energy E:

```
g(E) = Σₙₖ δ(E - εₙₖ)
```

Integrated DOS gives total electrons:

```
N = ∫_{-∞}^{E_F} g(E) dE
```

## Band Gap Estimation

From the DOS, we find the band gap by identifying:
- Valence band maximum (VBM): Highest occupied states
- Conduction band minimum (CBM): Lowest unoccupied states
- Gap = CBM - VBM

**Note**: PBE systematically underestimates band gaps. For accurate gaps, use HSE06 (see example 13).

## Output

```
Band gap (from DOS): 0.62 eV
Experimental band gap: 1.12 eV
Note: DFT-PBE underestimates band gaps
```

## Tips

1. **Dense K-Grid**: More k-points = smoother DOS
2. **Tetrahedron Method**: Most accurate for DOS
3. **Energy Range**: Adjust EMIN/EMAX if needed
4. **PDOS Analysis**: Use LORBIT=11 for orbital contributions

## Next Steps

Try `06_band_structure/` to plot band dispersion along high-symmetry paths.
