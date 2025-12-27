# 11 - Phonons

Calculate phonon properties of silicon using the finite displacement method.

## Key Concepts

- **Phonons**: Quantized lattice vibrations
- **Finite Displacement Method**: Calculate force constants from displaced configurations
- **Phonon DOS**: Distribution of vibrational frequencies
- **Thermal Properties**: Heat capacity, entropy, free energy from phonons

## Requirements

```bash
pip install phonopy
```

## VASP Parameters

### Force Calculation

```python
calc = Vasp(
    atoms=supercell,
    xc='PBE',
    encut=400,
    kpts=(4, 4, 4),  # Reduced for supercell
    ismear=0,
    sigma=0.05,
    ibrion=-1,       # No relaxation
    nsw=0,
)
```

### Key Settings

| Parameter | Value | Description |
|-----------|-------|-------------|
| `supercell` | 2×2×2 | Size for accurate dispersion |
| `displacement` | 0.01 Å | Small atomic displacement |
| `kpts` | Dense | Converged k-mesh for forces |

## Workflow

1. **Create supercell**: 2×2×2 or larger
2. **Generate displacements**: Using Phonopy
3. **Calculate forces**: VASP for each displaced structure
4. **Compute force constants**: From forces via Phonopy
5. **Analyze**: Dispersion, DOS, thermal properties

## Physical Background

### Force Constants

The dynamical matrix D(q) is constructed from force constants:

```
D_ij(q) = Σ_R Φ_ij(R) exp(iq·R) / √(M_i M_j)
```

where Φ_ij(R) are the force constants.

### Phonon Frequencies

Eigenvalues of D(q) give phonon frequencies:

```
D(q)|ε⟩ = ω²|ε⟩
```

### Thermal Properties

From phonon frequencies, calculate:

```
F_vib = kT Σ_q,j ln[2sinh(ℏω_qj/2kT)]

C_v = Σ_q,j kB (ℏω_qj/kT)² × [exp(ℏω/kT)/(exp(ℏω/kT)-1)²]
```

## Output

```
Phonon frequencies at Γ point (THz):
  Mode 1:   0.0000 THz (acoustic)
  Mode 2:   0.0000 THz (acoustic)
  Mode 3:   0.0000 THz (acoustic)
  Mode 4:  15.2300 THz (optical)
  Mode 5:  15.2300 THz (optical)
  Mode 6:  15.2300 THz (optical)

At T = 300 K:
  Heat capacity: 19.8 J/mol/K
  Experimental Cv: ~20 J/mol/K
```

## Tips

1. **Supercell Size**: Larger = more accurate dispersion
2. **Displacement**: 0.01 Å typical, use symmetry to reduce calculations
3. **K-points**: Well-converged forces are essential
4. **ENCUT**: Same as production calculations
5. **Symmetry**: Phonopy uses crystal symmetry to reduce displacements

## Common Issues

- **Imaginary frequencies**: Structure not fully relaxed
- **Acoustic modes**: Should be ~0 at Γ (check for small negative values)
- **Convergence**: Test supercell size and k-mesh

## Applications

- Zero-point energy corrections
- Thermodynamic properties
- Phase stability (free energy)
- Thermal conductivity (with phonon lifetimes)
- Raman/IR spectra

## Next Steps

Try `12_dft_plus_u/` for strongly correlated electron systems.
