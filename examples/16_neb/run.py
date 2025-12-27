#!/usr/bin/env python
"""
16 - Nudged Elastic Band (NEB)

This script demonstrates setting up and running NEB calculations to find
minimum energy paths and transition states for chemical reactions.

We use H diffusion on Cu(111) as a simple example.

Usage:
    python run.py

Note: NEB calculations are computationally expensive and require careful
      setup of initial and final states.
"""

from ase.build import fcc111, add_adsorbate
from ase.constraints import FixAtoms
from ase.neb import NEB, interpolate
import numpy as np
from vasp import Vasp

print("=" * 60)
print("Nudged Elastic Band (NEB) Calculations")
print("=" * 60)
print()

# =============================================================================
# Part 1: Create initial and final states
# =============================================================================

print("Part 1: Setting up endpoints")
print("-" * 40)
print()

# Create Cu(111) slab
slab = fcc111('Cu', size=(3, 3, 3), vacuum=10.0, a=3.6)

# Fix bottom layer
positions_z = slab.positions[:, 2]
z_sorted = np.sort(np.unique(np.round(positions_z, decimals=2)))
z_threshold = z_sorted[0] + 0.1
constraint = FixAtoms(mask=positions_z <= z_threshold)
slab.set_constraint(constraint)

# Initial state: H at fcc hollow site
initial = slab.copy()
add_adsorbate(initial, 'H', height=1.0, position='fcc')
initial.set_constraint(constraint)

# Final state: H at hcp hollow site
final = slab.copy()
add_adsorbate(final, 'H', height=1.0, position='hcp')
final.set_constraint(constraint)

print("System: H diffusion on Cu(111)")
print(f"  Slab: 3x3 Cu(111), {len(slab)} atoms")
print(f"  Initial: H at fcc hollow")
print(f"  Final: H at hcp hollow")
print()

# =============================================================================
# Part 2: Relax endpoints
# =============================================================================

print("Part 2: Relaxing endpoints")
print("-" * 40)
print()

# Relax initial state
print("Relaxing initial state (fcc site)...")
calc_initial = Vasp(
    label='neb/initial',
    atoms=initial,
    xc='PBE',
    encut=400,
    kpts=(4, 4, 1),
    ismear=1,
    sigma=0.1,
    ibrion=2,
    isif=2,
    nsw=50,
    ediffg=-0.02,
)

e_initial = calc_initial.potential_energy
print(f"  Initial state energy: {e_initial:.6f} eV")

# Relax final state
print("Relaxing final state (hcp site)...")
calc_final = Vasp(
    label='neb/final',
    atoms=final,
    xc='PBE',
    encut=400,
    kpts=(4, 4, 1),
    ismear=1,
    sigma=0.1,
    ibrion=2,
    isif=2,
    nsw=50,
    ediffg=-0.02,
)

e_final = calc_final.potential_energy
print(f"  Final state energy: {e_final:.6f} eV")
print()

# Get relaxed structures
initial_relaxed = calc_initial.atoms.copy()
final_relaxed = calc_final.atoms.copy()

print(f"  Reaction energy: {e_final - e_initial:.4f} eV")
print()

# =============================================================================
# Part 3: Create NEB images
# =============================================================================

print("Part 3: Creating NEB images")
print("-" * 40)
print()

n_images = 5  # Number of intermediate images

# Create list of images (initial + intermediates + final)
images = [initial_relaxed.copy()]
for i in range(n_images):
    images.append(initial_relaxed.copy())
images.append(final_relaxed.copy())

# Linear interpolation between endpoints
interpolate(images)

print(f"  Number of images: {n_images + 2} (including endpoints)")
print(f"  Interpolation: Linear")
print()

# =============================================================================
# Part 4: VASP NEB setup (conceptual)
# =============================================================================

print("Part 4: VASP NEB parameters")
print("-" * 40)
print()

print("For VASP's built-in NEB, set up the following directory structure:")
print()
print("  neb/")
print("  ├── INCAR       (NEB-specific parameters)")
print("  ├── KPOINTS")
print("  ├── POTCAR")
print("  ├── 00/POSCAR   (initial state)")
print("  ├── 01/POSCAR   (image 1)")
print("  ├── 02/POSCAR   (image 2)")
print("  ├── ...")
print("  └── 06/POSCAR   (final state)")
print()

neb_incar_params = """
# NEB-specific INCAR parameters:

# NEB method
IMAGES = 5        # Number of intermediate images
SPRING = -5.0     # Spring constant between images (eV/Å²)
LCLIMB = .TRUE.   # Climbing image for accurate TS

# Optimization
IBRION = 1        # Quasi-Newton (or 3 for damped MD)
POTIM = 0.1       # Step size
NSW = 100         # Max ionic steps
EDIFFG = -0.05    # Force convergence (eV/Å)

# Electronic
ISTART = 0
ICHARG = 2
ENCUT = 400
ISMEAR = 1
SIGMA = 0.1

# Parallelization
NCORE = 4         # Adjust for your system
"""

print("INCAR parameters for NEB:")
print(neb_incar_params)

# =============================================================================
# Part 5: Running with ASE NEB (alternative)
# =============================================================================

print("Part 5: ASE NEB alternative")
print("-" * 40)
print()

print("You can also run NEB using ASE's optimizer with VASP as calculator:")
print()

ase_neb_example = '''
from ase.neb import NEB
from ase.optimize import BFGS
from vasp import Vasp

# Create NEB object
neb = NEB(images, climb=True, parallel=False)

# Attach calculators to intermediate images
for image in images[1:-1]:
    calc = Vasp(
        atoms=image,
        xc='PBE',
        encut=400,
        kpts=(4, 4, 1),
        ismear=1,
        sigma=0.1,
    )
    image.calc = calc

# Optimize with BFGS
optimizer = BFGS(neb, trajectory='neb.traj')
optimizer.run(fmax=0.05)
'''

print(ase_neb_example)
print()

# =============================================================================
# Part 6: Analysis
# =============================================================================

print("Part 6: NEB analysis")
print("-" * 40)
print()

print("After NEB converges, analyze results:")
print()

analysis_code = '''
# Read energies along MEP
from ase.io import read

images = read('neb.traj', index=':')
energies = [img.get_potential_energy() for img in images]

# Find transition state
e_ts = max(energies)
ts_index = energies.index(e_ts)

# Calculate barrier
e_forward = e_ts - energies[0]    # Forward barrier
e_reverse = e_ts - energies[-1]   # Reverse barrier

print(f"Forward barrier: {e_forward:.3f} eV")
print(f"Reverse barrier: {e_reverse:.3f} eV")

# Verify TS (single imaginary frequency)
# Run frequency calculation on image with highest energy
'''

print(analysis_code)
print()

# =============================================================================
# Part 7: Expected results
# =============================================================================

print("Part 7: Expected results for H/Cu(111)")
print("-" * 40)
print()

print("H diffusion on Cu(111): fcc → hcp")
print()
print("  Literature values:")
print("    Forward barrier: ~0.15-0.20 eV")
print("    Reverse barrier: ~0.15-0.20 eV")
print("    (fcc and hcp sites are nearly isoenergetic)")
print()
print("  TS geometry:")
print("    H at bridge site between fcc and hcp hollows")
print()

# =============================================================================
# Part 8: Tips and troubleshooting
# =============================================================================

print("Part 8: Tips and troubleshooting")
print("-" * 40)
print()

print("Common issues and solutions:")
print()
print("  1. No saddle point found:")
print("     - Increase number of images")
print("     - Check initial path is reasonable")
print("     - Use IDPP or other improved interpolation")
print()
print("  2. Images bunching at endpoints:")
print("     - Adjust spring constant (SPRING)")
print("     - Try different interpolation scheme")
print()
print("  3. High barrier oscillates:")
print("     - Tighter force convergence")
print("     - More electronic steps")
print("     - Use climbing image (LCLIMB)")
print()
print("  4. TS verification:")
print("     - Run frequency calc on TS image")
print("     - Should have exactly one imaginary frequency")
print("     - Mode should point along reaction coordinate")
print()

# =============================================================================
# Part 9: Advanced NEB methods
# =============================================================================

print("Part 9: Advanced NEB methods")
print("-" * 40)
print()

print("VASP supports several NEB variants:")
print()
print("  Standard NEB:")
print("    IBRION = 1 or 3")
print("    IMAGES = N")
print("    SPRING = -5.0")
print()
print("  Climbing Image NEB (CI-NEB):")
print("    LCLIMB = .TRUE.")
print("    Better for accurate transition state")
print()
print("  VTST Tools extensions:")
print("    - Improved tangent estimates")
print("    - Dimer method for TS refinement")
print("    - Quasi-Newton optimization")
print("    See: http://theory.cm.utexas.edu/vtsttools/")
print()

# =============================================================================
# Summary
# =============================================================================

print("=" * 60)
print("Summary")
print("=" * 60)
print()
print("NEB workflow:")
print("  1. Relax initial and final states")
print("  2. Create interpolated images")
print("  3. Run NEB optimization")
print("  4. Analyze minimum energy path")
print("  5. Verify TS with frequency calculation")
print()
print("Key VASP parameters:")
print("  IMAGES = N      (number of intermediate images)")
print("  SPRING = -5.0   (spring constant)")
print("  LCLIMB = .TRUE. (climbing image for accurate TS)")
print("  EDIFFG = -0.05  (force convergence)")
print()
print("Applications:")
print("  - Diffusion barriers")
print("  - Reaction mechanisms")
print("  - Defect migration")
print("  - Catalytic pathways")
print()
print("This completes the NEB tutorial!")
