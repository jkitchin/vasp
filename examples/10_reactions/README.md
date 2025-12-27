# 10 - Reactions

Calculate reaction energetics: dissociation energies, reaction enthalpies, and transition state concepts.

## Key Concepts

- **Reaction Energy**: Energy difference between products and reactants
- **Dissociation Energy**: Energy to break a bond
- **Activation Energy**: Barrier height (transition state - reactants)
- **Potential Energy Surface (PES)**: Energy as function of coordinates

## VASP Parameters

### Molecule Calculations

```python
h2 = molecule('H2')
h2.center(vacuum=8.0)

calc = Vasp(
    atoms=h2,
    xc='PBE',
    encut=400,
    kpts=(1, 1, 1),
    ismear=0,
    sigma=0.05,
    ispin=2,
)
```

### Triplet O2 (Important!)

```python
calc = Vasp(
    atoms=o2,
    xc='PBE',
    encut=400,
    kpts=(1, 1, 1),
    ispin=2,
    magmom=[1.0, 1.0],  # Triplet state
)
```

## Physical Background

### Reaction Energy

```
ΔE = Σ E(products) - Σ E(reactants)
```

- **ΔE < 0**: Exothermic (favorable)
- **ΔE > 0**: Endothermic (unfavorable)

### Example: Combustion

```
2H2 + O2 → 2H2O
ΔE = 2×E(H2O) - 2×E(H2) - E(O2) ≈ -5.0 eV
```

### Activation Energy

```
Ea = E(transition state) - E(reactants)
```

The Arrhenius equation relates this to rate:
```
k = A × exp(-Ea / kT)
```

## Transition State Methods

### NEB (Nudged Elastic Band)

```python
# VASP NEB settings
IMAGES = 5          # Number of intermediate images
SPRING = -5.0       # Spring constant
IBRION = 1          # Optimizer
POTIM = 0.1         # Step size
LCLIMB = True       # Climbing image for accurate TS
```

### Workflow

1. Optimize reactant and product structures
2. Create initial path (linear interpolation)
3. Run NEB to find MEP (minimum energy path)
4. Identify transition state (highest point)
5. Verify with frequency calculation (one imaginary frequency)

## Output

```
Key reaction energetics:
  H2 dissociation: 4.48 eV (exp: 4.52 eV)
  H2 combustion: -2.48 eV/H2O (exp: ~-2.5 eV)
```

## Tips

1. **Spin States**: Use correct spin multiplicity (O2 is triplet)
2. **Box Size**: Large enough vacuum for molecules
3. **Zero-Point Energy**: Include for accurate energetics
4. **Entropy**: Important for gas-phase reactions
5. **NEB Images**: 5-9 images typically sufficient

## Common Corrections

### Zero-Point Energy (ZPE)

```
ZPE = (1/2) × Σ ℏω
```

Include for accurate ΔE and Ea.

### Thermodynamic Corrections

For T > 0:
```
G = E + ZPE + H_vib - T×S
```

## Applications

- Catalytic reaction mechanisms
- Bond dissociation energies
- Reaction selectivity
- Microkinetic modeling

## Next Steps

Try `11_phonons/` for vibrational properties and zero-point energies.
