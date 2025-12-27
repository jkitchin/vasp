# 09 - Adsorption

Calculate adsorption energies for CO on Pt(111) at different sites.

## Key Concepts

- **Adsorption Energy**: Energy gained upon binding
- **Adsorption Sites**: Ontop, bridge, hollow (fcc/hcp)
- **Coverage**: Number of adsorbates per surface site
- **Reference State**: Isolated molecule and clean surface

## VASP Parameters

### Gas-Phase Molecule

```python
co = molecule('CO')
co.center(vacuum=8.0)

calc = Vasp(
    atoms=co,
    xc='PBE',
    encut=400,
    kpts=(1, 1, 1),   # Gamma only
    ismear=0,
    sigma=0.05,
    ispin=2,          # Open shell
)
```

### Surface + Adsorbate

```python
from ase.build import fcc111, add_adsorbate

slab = fcc111('Pt', size=(2, 2, 4), vacuum=12.0)
add_adsorbate(slab, 'CO', height=1.9, position='ontop')

calc = Vasp(
    atoms=slab,
    xc='PBE',
    encut=400,
    kpts=(4, 4, 1),
    ismear=1,
    sigma=0.1,
    ldipol=True,
    idipol=3,
    isif=2,
    ibrion=2,
    nsw=50,
    ediffg=-0.03,
)
```

## Adsorption Sites (FCC 111)

| Site | Description | Coordination |
|------|-------------|--------------|
| Ontop | Above single atom | 1 |
| Bridge | Between two atoms | 2 |
| FCC | 3-fold hollow (no atom below) | 3 |
| HCP | 3-fold hollow (atom below) | 3 |

## Physical Background

### Adsorption Energy

```
E_ads = E(surface+adsorbate) - E(clean surface) - E(gas molecule)
```

- **Negative E_ads**: Favorable (exothermic) adsorption
- **Positive E_ads**: Unfavorable (endothermic)

### Coverage Definition

Coverage θ is fraction of surface sites occupied:
- θ = 1 ML: One adsorbate per surface atom
- θ = 0.25 ML: One adsorbate per 2×2 cell

### CO on Pt(111)

Experimental findings:
- Preferred site: Atop at low coverage
- Adsorption energy: ~-1.4 to -1.5 eV
- C-O stretch: ~2080 cm⁻¹
- Pt-CO distance: ~1.85 Å

## Output

```
Summary: CO/Pt(111) Adsorption

Site       E_ads (eV)   Stability
ontop      -1.52        Most stable
bridge     -1.38        ΔE = 0.14 eV
fcc        -1.35        ΔE = 0.17 eV
hcp        -1.33        ΔE = 0.19 eV
```

Note: DFT-PBE often overestimates CO binding. Consider BEEF-vdW or RPBE for more accurate energetics.

## Tips

1. **Reference Calculations**: Use same parameters for all
2. **Adsorbate Orientation**: Check C-down vs O-down for CO
3. **Site Search**: Try all high-symmetry sites
4. **Relaxation**: Allow adsorbate and top layers to relax
5. **Spin Polarization**: Check if molecule has unpaired electrons

## Common Issues

- **CO puzzle**: Standard DFT predicts wrong site preference for CO on some metals
- **Van der Waals**: Important for physisorbed species
- **Zero-point energy**: Can shift relative stabilities

## Applications

- Heterogeneous catalysis
- Surface poisoning
- Electrochemistry
- Sensor design

## Next Steps

Try `10_reactions/` for reaction energy barriers.
