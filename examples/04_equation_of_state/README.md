# 04 - Equation of State

Calculate the equation of state (E vs V) and determine the bulk modulus of silicon.

## Key Concepts

- **Equation of State**: Relationship between energy and volume
- **Bulk Modulus**: Resistance to uniform compression
- **Birch-Murnaghan EOS**: Standard fitting function for solids

## VASP Parameters

This example uses static calculations at different volumes:

```python
calc = Vasp(
    atoms=atoms,
    xc='PBE',
    encut=400,
    kpts=(8, 8, 8),
    ismear=1,
    sigma=0.1,
)
```

## Physical Background

The bulk modulus B₀ is defined as:

```
B₀ = -V (∂P/∂V)_T = V (∂²E/∂V²)
```

We fit E(V) data to the Birch-Murnaghan equation:

```
E(V) = E₀ + (9V₀B₀/16) × {[(V₀/V)^(2/3) - 1]³B₀' + [(V₀/V)^(2/3) - 1]² × [6 - 4(V₀/V)^(2/3)]}
```

where:
- E₀: Equilibrium energy
- V₀: Equilibrium volume
- B₀: Bulk modulus
- B₀': Pressure derivative of bulk modulus

## Output

```
Fitted parameters:
  Equilibrium volume: 40.04 Å³
  Equilibrium energy: -10.845 eV
  Bulk modulus: 88.5 GPa
  B' (pressure derivative): 4.2

Experimental bulk modulus: 99 GPa
```

## Tips

1. **Volume Range**: Sample ±10% around equilibrium
2. **Number of Points**: 7-11 points usually sufficient
3. **Dense K-points**: Use converged k-grid for each volume
4. **Consistent Settings**: Same ENCUT and k-points for all volumes

## Applications

- Determine equilibrium lattice constant
- Calculate elastic constants
- Study pressure effects
- Validate pseudopotentials

## Next Steps

Try `05_density_of_states/` to learn about electronic structure analysis.
