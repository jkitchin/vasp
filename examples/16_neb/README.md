# 16 - Nudged Elastic Band (NEB)

Find minimum energy paths and transition states using the NEB method.

## Key Concepts

- **NEB**: Chain of images connected by springs along reaction path
- **MEP**: Minimum Energy Path from reactants to products
- **Transition State**: Saddle point on the potential energy surface
- **Climbing Image**: Improved method for accurate TS location

## VASP Parameters

### NEB Setup

```
# INCAR for NEB calculation
IMAGES = 5        # Intermediate images
SPRING = -5.0     # Spring constant (eV/Å²)
LCLIMB = .TRUE.   # Climbing image NEB

IBRION = 1        # Quasi-Newton optimizer
POTIM = 0.1       # Step size
NSW = 100         # Max steps
EDIFFG = -0.05    # Force convergence (eV/Å)
```

### Directory Structure

```
neb/
├── INCAR
├── KPOINTS
├── POTCAR
├── 00/POSCAR   # Initial state
├── 01/POSCAR   # Image 1
├── 02/POSCAR   # Image 2
├── ...
└── 06/POSCAR   # Final state
```

## Workflow

1. **Prepare endpoints**
   - Relax initial state (reactants)
   - Relax final state (products)

2. **Create images**
   - Linear interpolation (simple)
   - IDPP interpolation (better for complex paths)

3. **Run NEB**
   - VASP with IMAGES tag
   - Or ASE NEB with VASP calculator

4. **Analyze**
   - Extract energies along path
   - Identify transition state
   - Calculate barriers

5. **Verify**
   - Frequency calculation at TS
   - Single imaginary frequency expected

## Physical Background

### Spring Force

NEB adds spring forces to keep images evenly spaced:

```
F_spring = k × (|R_i+1 - R_i| - |R_i - R_i-1|) × τ
```

where τ is the tangent to the path.

### Climbing Image

The highest-energy image "climbs" toward the saddle point:

```
F_CI = F_true - 2(F_true · τ)τ
```

This gives the exact transition state.

## Example: H Diffusion on Cu(111)

```
Initial: H at fcc hollow
Final: H at hcp hollow
Path: fcc → bridge (TS) → hcp

Forward barrier: ~0.15-0.20 eV
```

## Tips

1. **Number of Images**: 5-9 typically sufficient
2. **Spring Constant**: -5 to -10 eV/Å² works for most cases
3. **Convergence**: Use tighter EDIFFG for accurate barriers
4. **Initial Path**: Quality of interpolation affects convergence

## Common Issues

| Problem | Solution |
|---------|----------|
| Images bunch at endpoints | Adjust spring constant |
| No saddle point | More images, better interpolation |
| Slow convergence | Try different optimizer (IBRION) |
| Wrong TS | Verify with frequency calculation |

## TS Verification

Run frequency calculation on TS image:

```python
calc = Vasp(
    atoms=ts_image,
    ibrion=5,       # Finite differences
    nfree=2,        # Central differences
    potim=0.015,    # Displacement
    nsw=1,
)
```

Check for:
- Exactly one imaginary frequency
- Mode pointing along reaction coordinate

## Applications

- Surface diffusion barriers
- Heterogeneous catalysis mechanisms
- Defect migration in solids
- Molecular conformational changes
- Reaction selectivity

## Advanced Methods

### VTST Tools

Enhanced NEB implementation:
- Improved tangent estimates
- Dimer method for TS refinement
- Better optimization schemes

Website: http://theory.cm.utexas.edu/vtsttools/

### Dimer Method

Alternative for finding TS without knowing products:
- Start near suspected TS
- Converge to nearest saddle point
- Doesn't require final state
