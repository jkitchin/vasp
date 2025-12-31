# Beginner Tutorials

These tutorials cover the fundamentals of DFT calculations with VASP.

## Overview

| Tutorial | Topic | Description | Notebook |
|----------|-------|-------------|----------|
| 01 | Getting Started | First VASP calculation | [View](https://github.com/jkitchin/vasp/blob/master/docs/tutorials/01_getting_started.ipynb) |
| 02 | Convergence | ENCUT and k-point convergence | [View](https://github.com/jkitchin/vasp/blob/master/docs/tutorials/02_convergence.ipynb) |
| 03 | Relaxation | Structure optimization | [View](https://github.com/jkitchin/vasp/blob/master/docs/tutorials/03_relaxation.ipynb) |
| 04 | Equation of State | Bulk modulus calculation | [View](https://github.com/jkitchin/vasp/blob/master/docs/tutorials/04_equation_of_state.ipynb) |
| 05 | Density of States | Electronic structure | [View](https://github.com/jkitchin/vasp/blob/master/docs/tutorials/05_density_of_states.ipynb) |

## Prerequisites

- VASP installed and configured
- Basic Python knowledge
- Familiarity with atomic simulation concepts

## Topics Covered

### 01 - Getting Started

Learn how to:
- Create an atomic structure with ASE
- Set up a VASP calculator
- Run a single-point energy calculation
- Read results

### 02 - Convergence Testing

Learn how to:
- Test ENCUT convergence
- Test k-point convergence
- Determine production settings

### 03 - Structure Relaxation

Learn how to:
- Relax atomic positions (ISIF=2)
- Relax cell and positions (ISIF=3)
- Set force convergence criteria

### 04 - Equation of State

Learn how to:
- Calculate energy vs. volume
- Fit Birch-Murnaghan equation
- Determine bulk modulus

### 05 - Density of States

Learn how to:
- Perform SCF + non-SCF workflow
- Calculate total and projected DOS
- Estimate band gap from DOS
