# 06 - Band Structure

Calculate the electronic band structure of silicon along high-symmetry paths.

## Key Concepts

- **Band Structure**: Energy dispersion E(k) along k-paths
- **High-Symmetry Points**: Special k-points in the Brillouin zone
- **Direct vs Indirect Gap**: Where VBM and CBM occur in k-space
- **Two-Step Calculation**: SCF → non-SCF bands

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
    lwave=True,
    lcharg=True,
)
```

### Step 2: Non-Self-Consistent Bands

```python
calc_bands = Vasp(
    atoms=atoms,
    xc='PBE',
    encut=400,
    kpts=kpath,           # K-points along path
    reciprocal=True,       # Fractional coordinates
    ismear=0,              # Gaussian smearing
    sigma=0.05,
    icharg=11,             # Read charge from file
    lorbit=11,             # Orbital character
)
```

## High-Symmetry Points (FCC)

| Point | Coordinates | Description |
|-------|-------------|-------------|
| Γ | (0, 0, 0) | Zone center |
| X | (0.5, 0, 0.5) | Zone face |
| L | (0.5, 0.5, 0.5) | Zone corner |
| K | (0.375, 0.375, 0.75) | Zone edge |
| U | (0.625, 0.25, 0.625) | Zone edge |

Common path: L → Γ → X → U|K → Γ

## Output Files

- **EIGENVAL**: Eigenvalues at each k-point
- **PROCAR**: Orbital projections (with LORBIT=11)

## Physical Background

Silicon is an indirect gap semiconductor:
- **VBM**: Located at Γ point
- **CBM**: Located near X point (actually along Γ-X)
- **Indirect gap**: ~1.12 eV (experimental)
- **Direct gap at Γ**: ~3.4 eV

DFT-PBE underestimates the gap by ~0.5 eV.

## Output

```
Valence band maximum: -0.02 eV (relative to Fermi)
Conduction band minimum: 0.60 eV
Band gap: 0.62 eV
Experimental band gap: 1.12 eV
```

## Tips

1. **K-Path Density**: Use ~30-50 points per segment
2. **Gaussian Smearing**: ISMEAR=0 for band calculations
3. **Charge Reading**: Always copy CHGCAR from SCF
4. **Orbital Character**: LORBIT=11 for fat bands

## Applications

- Identify semiconductor type (direct/indirect)
- Analyze orbital contributions to bands
- Study band dispersion and effective mass
- Compare different exchange-correlation functionals

## Next Steps

Try `07_magnetism/` for spin-polarized calculations.
