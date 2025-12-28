# 17 - Molecular Vibrations

Calculate vibrational frequencies and modes for molecules using finite differences.

## Key Concepts

- **Vibrational Frequencies**: Normal modes of molecular vibration
- **Finite Differences**: Calculate Hessian from displaced configurations
- **IR/Raman Activity**: Dipole and polarizability changes with vibration
- **Zero-Point Energy**: Quantum vibrational ground state energy

## VASP Parameters

### Force Calculation

```python
calc = Vasp(
    atoms=molecule,
    xc='PBE',
    encut=400,
    kpts=(1, 1, 1),  # Gamma point for molecules
    ismear=0,
    sigma=0.01,
    ibrion=5,        # Finite differences for vibrations
    nfree=2,         # Central differences (recommended)
    potim=0.015,     # Displacement in Angstrom
)
```

### Key Settings

| Parameter | Value | Description |
|-----------|-------|-------------|
| `ibrion` | 5 | Finite differences for Hessian |
| `nfree` | 2 | Central differences (more accurate) |
| `potim` | 0.015 | Step size for displacements (Å) |
| `ismear` | 0 | Gaussian smearing for molecules |

### Alternative: DFPT Method

```python
calc = Vasp(
    atoms=molecule,
    ibrion=7,        # DFPT (density functional perturbation theory)
    # or ibrion=8    # DFPT including symmetry
)
```

## Physical Background

### Harmonic Approximation

Potential energy surface expanded to second order:

```
E(R) ≈ E₀ + ½ ΣᵢⱼΣₐᵦ Φᵢⱼ,ₐᵦ uᵢₐ uⱼᵦ
```

where Φ is the Hessian (force constant matrix).

### Normal Modes

Diagonalization of mass-weighted Hessian:

```
D = M⁻¹/² Φ M⁻¹/²
D |ε⟩ = ω² |ε⟩
```

gives frequencies ω and mode vectors ε.

### Zero-Point Energy

```
E_ZPE = ½ Σₙ ℏωₙ
```

## Output

```
Vibrational Frequencies for H2O:
  Mode 1: 1595.3 cm⁻¹  (bending)
  Mode 2: 3657.1 cm⁻¹  (symmetric stretch)
  Mode 3: 3756.0 cm⁻¹  (asymmetric stretch)

Experimental values:
  Bending:    1595 cm⁻¹
  Sym. str.:  3657 cm⁻¹
  Asym. str.: 3756 cm⁻¹

Zero-point energy: 0.558 eV (55.8 kJ/mol)
```

## Tips

1. **Box Size**: Large vacuum (>10 Å) for isolated molecules
2. **Step Size**: 0.015 Å typically good; test for convergence
3. **Symmetry**: Use ISYM=0 for accurate frequencies
4. **Convergence**: Tight EDIFF (1e-7 or better) for accurate forces
5. **Imaginary Frequencies**: Indicate saddle point (transition state)

## Common Issues

- **Translation/Rotation Modes**: Should be ~0 cm⁻¹ (6 modes for nonlinear)
- **Large Box**: Use 15+ Å vacuum for accurate dipole derivatives
- **POTIM**: Too large causes anharmonic errors; too small causes numerical noise

## Applications

- IR/Raman spectra prediction
- Thermochemistry (ZPE corrections)
- Isotope effects
- Transition state verification (one imaginary frequency)

## Next Steps

Try `18_3d_visualization/` for visualizing charge density, potential, and ELF.
