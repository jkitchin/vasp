# VASP Tutorials

Tutorials are now located in [`docs/tutorials/`](../docs/tutorials/).

## Running Tutorials

```bash
# Run all tutorials
make tutorials

# Run a single tutorial by number or name
make tutorial T=01
make tutorial T=getting_started

# Clear notebook outputs
make tutorials-clean
```

## Tutorial List

| Tutorial | Topic | Description |
|----------|-------|-------------|
| `01_getting_started` | First calculation | Energy of a silicon crystal |
| `02_convergence` | Convergence testing | ENCUT and k-point convergence |
| `03_relaxation` | Structure optimization | Relaxing atomic positions and cell |
| `04_equation_of_state` | Bulk modulus | Fitting energy vs volume |
| `05_density_of_states` | Electronic structure | Total and projected DOS |
| `06_band_structure` | Band structure | Band structure along high-symmetry path |
| `07_magnetism` | Spin polarization | Magnetic moments in Fe and NiO |
| `08_surfaces` | Slab calculations | Surface energy of Cu(111) |
| `09_adsorption` | Molecules on surfaces | CO on Pt(111) |
| `10_reactions` | Reaction energetics | Dissociation energy and barriers |
| `11_phonons` | Vibrational properties | Phonon dispersion with Phonopy |
| `12_dft_plus_u` | Strongly correlated | DFT+U for FeO |
| `13_hybrid_functionals` | HSE06 | Band gap with hybrid functional |
| `14_van_der_waals` | Dispersion corrections | Graphite interlayer binding |
| `15_workflows` | Automated pipelines | Using recipes for high-throughput |
| `16_neb` | Nudged elastic band | Transition states and reaction barriers |
| `17_vibrations` | Molecular vibrations | IR/Raman frequencies for molecules |
| `18_3d_visualization` | Volumetric data | Charge density, potential, ELF |
| `19_pseudopotentials` | POTCAR selection | Choosing the right pseudopotentials |
| `20_interactive_mode` | Interactive VASP | Persistent process for fast optimization |
| `21_cluster_expansion` | Cluster expansion | Ni-Al CE with ICET and SQS |

See the full documentation at https://kitchingroup.cheme.cmu.edu/vasp/
