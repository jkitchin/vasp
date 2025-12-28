# Example 21: Cluster Expansion with ICET

## Overview

This example demonstrates using ICET (Integrated Cluster Expansion Toolkit) to build a cluster expansion for the Ni-Al system. We use a minimal training set to illustrate the workflow.

## What is a Cluster Expansion?

A cluster expansion represents the energy of an alloy as a sum over clusters:

```
E(σ) = Σ_α m_α J_α Π_α(σ)
```

Where:
- `σ` is the atomic configuration
- `J_α` are effective cluster interactions (ECIs)
- `Π_α(σ)` are cluster correlation functions

## Requirements

```bash
pip install icet
```

## Workflow

1. **Define parent lattice** - FCC structure with Ni/Al sites
2. **Generate training structures** - Pure elements + ordered compounds
3. **Calculate DFT energies** - Using VASP
4. **Fit cluster expansion** - Determine ECIs
5. **Generate SQS** - Special Quasirandom Structure for random alloy

## Training Set

We use a minimal set of 5 structures:
- Pure Ni (FCC)
- Pure Al (FCC)
- NiAl (B2, CsCl-type)
- Ni₃Al (L1₂)
- NiAl₃ (L1₂)

## Files

- `tutorial.ipynb` - Interactive notebook
- `README.md` - This documentation

## References

- [ICET Documentation](https://icet.materialsmodeling.org/)
- [Cluster Expansion Theory](https://doi.org/10.1103/RevModPhys.74.11)
