#!/usr/bin/env python
"""
15 - Workflows

This script demonstrates using the vasp-ase recipe system for automated
calculation workflows, including multi-step relaxations and complex
pipelines.

Usage:
    python run.py
"""

from ase.build import bulk, fcc111
import numpy as np
from vasp import Vasp
from vasp.runners import MockRunner, MockResults
from vasp.recipes.core import static_job, relax_job, double_relax_flow, VaspResult
from vasp.recipes.slabs import slab_static_job, slab_relax_job, bulk_to_slabs_flow
from vasp.recipes.decorators import job, flow, subflow

print("=" * 60)
print("Automated Workflows with Recipes")
print("=" * 60)
print()

# =============================================================================
# Part 1: Recipe basics
# =============================================================================

print("Part 1: Recipe system overview")
print("-" * 40)
print()

print("The recipe system provides reusable calculation patterns:")
print()
print("  @job: Single VASP calculation")
print("    - static_job: Single-point energy")
print("    - relax_job: Geometry optimization")
print()
print("  @flow: Multi-step workflow")
print("    - double_relax_flow: Coarse → fine relaxation")
print("    - bulk_to_slabs_flow: Bulk → slab generation → calculation")
print()
print("  @subflow: Parallel sub-calculations")
print("    - Process multiple structures in parallel")
print()

# =============================================================================
# Part 2: Using MockRunner for demonstration
# =============================================================================

print("Part 2: Demonstration with MockRunner")
print("-" * 40)
print()

print("Using MockRunner to demonstrate workflow without VASP:")
print()

# Create mock runner that returns predefined results
mock_results = MockResults(
    energy=-10.5,
    forces=np.zeros((2, 3)),
    stress=np.zeros(6),
)
runner = MockRunner(results=mock_results)

# Create silicon structure
si = bulk('Si', 'diamond', a=5.43)

# Run static calculation
print("Running static_job...")
result = static_job(si, runner=runner)
print(f"  Energy: {result.energy:.4f} eV")
print(f"  Converged: {result.converged}")
print()

# =============================================================================
# Part 3: Relaxation workflow
# =============================================================================

print("Part 3: Relaxation workflow")
print("-" * 40)
print()

print("Running relax_job...")
relax_result = relax_job(si, runner=runner, relax_cell=True)
print(f"  Energy: {relax_result.energy:.4f} eV")
print(f"  Parameters used: {list(relax_result.parameters.keys())[:5]}...")
print()

print("Running double_relax_flow (coarse → fine)...")
double_result = double_relax_flow(si, runner=runner)
print(f"  Final energy: {double_result.energy:.4f} eV")
print()

# =============================================================================
# Part 4: Slab workflow
# =============================================================================

print("Part 4: Slab workflow")
print("-" * 40)
print()

# Create copper bulk
cu = bulk('Cu', 'fcc', a=3.6)

print("Running bulk_to_slabs_flow for Cu...")
print("  This workflow:")
print("    1. Relaxes bulk structure")
print("    2. Generates slabs for specified Miller indices")
print("    3. Relaxes each slab")
print()

slab_results = bulk_to_slabs_flow(
    cu,
    runner=runner,
    miller_indices=[(1, 1, 1)],
)

print(f"  Bulk energy: {slab_results['bulk'].energy:.4f} eV")
print(f"  Number of slabs: {len(slab_results['slabs'])}")
for i, slab_result in enumerate(slab_results['slabs']):
    print(f"  Slab {i+1} energy: {slab_result.energy:.4f} eV")
print()

# =============================================================================
# Part 5: Custom recipe
# =============================================================================

print("Part 5: Creating custom recipes")
print("-" * 40)
print()

@job
def convergence_test_job(atoms, runner=None, encuts=None, **kwargs):
    """Test ENCUT convergence."""
    if encuts is None:
        encuts = [300, 400, 500, 600]

    results = []
    for encut in encuts:
        calc = Vasp(
            atoms=atoms,
            runner=runner,
            encut=encut,
            **kwargs
        )
        energy = calc.potential_energy
        results.append({'encut': encut, 'energy': energy})

    return results

print("Custom recipe: convergence_test_job")
print()

# Mock different energies for each ENCUT
test_results = convergence_test_job(
    si,
    runner=runner,
    encuts=[300, 400, 500],
    kpts=(4, 4, 4),
)
print("  Convergence test results:")
for r in test_results:
    print(f"    ENCUT={r['encut']}: E={r['energy']:.4f} eV")
print()

# =============================================================================
# Part 6: Workflow composition
# =============================================================================

print("Part 6: Composing workflows")
print("-" * 40)
print()

@flow
def full_analysis_flow(atoms, runner=None, **kwargs):
    """Complete analysis: relax → DOS → bands."""
    print("  Step 1: Relaxation")
    relax = relax_job(atoms, runner=runner, **kwargs)

    print("  Step 2: Static for DOS")
    static = static_job(relax.atoms, runner=runner, **kwargs)

    print("  Step 3: Band structure (would be added here)")
    # bands = band_structure_job(relax.atoms, runner=runner, **kwargs)

    return {
        'relaxed': relax,
        'static': static,
        # 'bands': bands,
    }

print("Custom flow: full_analysis_flow")
print()
results = full_analysis_flow(si, runner=runner)
print()
print(f"  Relaxed energy: {results['relaxed'].energy:.4f} eV")
print(f"  Static energy: {results['static'].energy:.4f} eV")
print()

# =============================================================================
# Part 7: Workflow engine integration
# =============================================================================

print("Part 7: Workflow engine integration")
print("-" * 40)
print()

print("The recipe decorators support workflow engines:")
print()
print("  Environment variable: VASP_WORKFLOW_ENGINE")
print()
print("  Supported engines:")
print("    - prefect: Prefect workflows")
print("    - dask: Dask distributed")
print("    - parsl: Parsl parallel")
print("    - covalent: Covalent cloud")
print("    - jobflow: Materials Project jobflow")
print()
print("  When an engine is set:")
print("    - @job becomes the engine's task decorator")
print("    - @flow becomes the engine's flow decorator")
print("    - Enables parallel execution, retries, logging")
print()

print("Example with Prefect:")
print('''
  export VASP_WORKFLOW_ENGINE=prefect

  from vasp.recipes import static_job, relax_job

  # Now decorated with @prefect.task / @prefect.flow
  # Can use Prefect features:
  # - Automatic retries
  # - Task caching
  # - Distributed execution
  # - Web UI monitoring
''')
print()

# =============================================================================
# Part 8: Real-world example structure
# =============================================================================

print("Part 8: Example production workflow")
print("-" * 40)
print()

example_workflow = '''
from vasp.recipes import relax_job, static_job, double_relax_flow
from vasp.recipes.phonons import phonon_flow
from vasp.runners import SlurmRunner

# Configure for HPC
runner = SlurmRunner(
    partition='compute',
    nodes=2,
    time='4:00:00',
)

# Material discovery workflow
def screen_material(atoms, name):
    """Screen a material for stability and properties."""

    # 1. Relax structure
    print(f"Relaxing {name}...")
    relaxed = double_relax_flow(atoms, runner=runner)

    if relaxed.energy is None:
        return None

    # 2. Check dynamical stability
    print(f"Checking phonons for {name}...")
    phonons = phonon_flow(relaxed.atoms, runner=runner)

    if any(f < -0.5 for f in phonons.frequencies):
        print(f"{name}: Dynamically unstable!")
        return None

    # 3. Calculate electronic properties
    print(f"Electronic properties for {name}...")
    electronic = static_job(
        relaxed.atoms,
        runner=runner,
        lorbit=11,
        nedos=2000,
    )

    return {
        'name': name,
        'structure': relaxed.atoms,
        'energy': relaxed.energy,
        'band_gap': electronic.band_gap,
        'phonons': phonons,
    }

# Screen multiple materials
materials = [bulk('Si'), bulk('Ge'), bulk('GaAs')]
names = ['Si', 'Ge', 'GaAs']

results = []
for atoms, name in zip(materials, names):
    result = screen_material(atoms, name)
    if result:
        results.append(result)
'''

print("Production workflow example:")
print(example_workflow)
print()

# =============================================================================
# Summary
# =============================================================================

print("=" * 60)
print("Summary")
print("=" * 60)
print()

print("Recipe system features:")
print("  - Reusable calculation patterns")
print("  - Consistent parameter handling")
print("  - Result dataclasses for type safety")
print("  - Workflow engine integration")
print("  - Easy composition of complex workflows")
print()

print("Available recipes:")
print("  Core: static_job, relax_job, double_relax_flow")
print("  Slabs: slab_static_job, slab_relax_job, bulk_to_slabs_flow")
print("  Phonons: phonon_job, phonon_flow (requires phonopy)")
print()

print("Best practices:")
print("  - Use recipes for reproducibility")
print("  - Set workflow engine for HPC")
print("  - Create custom recipes for repeated tasks")
print("  - Use MockRunner for testing")
print()

print("Congratulations! You've completed all examples.")
print("For more information, see the documentation.")
