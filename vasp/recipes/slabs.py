"""Slab and surface calculation recipes.

Provides workflows for surface calculations including:
- Slab generation from bulk
- Surface relaxation
- Adsorption calculations
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from .core import VaspResult, relax_job, static_job
from .decorators import flow, job, subflow

if TYPE_CHECKING:
    from ase import Atoms

    from ..runners.base import Runner


@dataclass
class SlabResult(VaspResult):
    """Result container for slab calculations.

    Additional attributes:
        miller_indices: Miller indices of the surface.
        surface_energy: Surface energy in J/m².
        work_function: Work function in eV.
        layers: Number of atomic layers.
    """
    miller_indices: tuple[int, int, int] | None = None
    surface_energy: float | None = None
    work_function: float | None = None
    layers: int = 0


def make_slabs_from_bulk(
    bulk_atoms: Atoms,
    miller_indices: list[tuple[int, int, int]] | None = None,
    min_slab_size: float = 10.0,
    min_vacuum_size: float = 15.0,
    max_normal_search: int = 1,
    center_slab: bool = True,
    symmetrize: bool = False,
) -> list[Atoms]:
    """Generate slabs from bulk structure.

    Args:
        bulk_atoms: Bulk atomic structure.
        miller_indices: List of Miller indices. If None, uses common surfaces.
        min_slab_size: Minimum slab thickness in Å.
        min_vacuum_size: Vacuum thickness in Å.
        max_normal_search: Max search depth for surface normals.
        center_slab: Center slab in cell.
        symmetrize: Make symmetric slabs.

    Returns:
        List of slab Atoms objects.
    """
    from ase.build import surface

    if miller_indices is None:
        # Common low-index surfaces
        miller_indices = [
            (1, 0, 0),
            (1, 1, 0),
            (1, 1, 1),
        ]

    slabs = []
    for indices in miller_indices:
        # Estimate layers needed for min_slab_size
        # This is approximate; actual thickness depends on structure
        layers = max(4, int(min_slab_size / 2.5))

        slab = surface(
            bulk_atoms,
            indices,
            layers=layers,
            vacuum=min_vacuum_size / 2,  # vacuum on each side
            periodic=True,
        )

        if center_slab:
            # Center in z direction
            slab.center(vacuum=min_vacuum_size / 2, axis=2)

        # Store miller indices in info
        slab.info['miller_indices'] = indices
        slab.info['layers'] = layers

        slabs.append(slab)

    return slabs


@job
def slab_static_job(
    slab: Atoms,
    runner: Runner | None = None,
    dipole_correction: bool = True,
    **calc_kwargs
) -> SlabResult:
    """Static calculation for a slab.

    Includes dipole correction by default for accurate surface properties.

    Args:
        slab: Slab atomic structure.
        runner: Job runner.
        dipole_correction: Apply dipole correction (LDIPOL).
        **calc_kwargs: Additional VASP parameters.

    Returns:
        SlabResult with energies.
    """
    from ..calculator import Vasp
    from ..parameters import get_preset

    params = get_preset('static')

    if dipole_correction:
        params['ldipol'] = True
        params['idipol'] = 3  # Correct in z direction

    params.update(calc_kwargs)

    calc = Vasp(
        atoms=slab,
        runner=runner,
        **params
    )

    slab.calc = calc
    energy = slab.get_potential_energy()

    return SlabResult(
        atoms=slab.copy(),
        energy=energy,
        forces=calc.results.get('forces'),
        stress=calc.results.get('stress'),
        fermi_level=calc.results.get('fermi_level'),
        parameters=calc.parameters.copy(),
        directory=calc.directory,
        miller_indices=slab.info.get('miller_indices'),
        layers=slab.info.get('layers', 0),
    )


@job
def slab_relax_job(
    slab: Atoms,
    runner: Runner | None = None,
    fix_bottom_layers: int = 2,
    relax_cell: bool = False,
    fmax: float = 0.02,
    dipole_correction: bool = True,
    **calc_kwargs
) -> SlabResult:
    """Relaxation for a slab with fixed bottom layers.

    Args:
        slab: Slab atomic structure.
        runner: Job runner.
        fix_bottom_layers: Number of bottom layers to fix.
        relax_cell: If True, relax in-plane cell parameters.
        fmax: Force convergence in eV/Å.
        dipole_correction: Apply dipole correction.
        **calc_kwargs: Additional VASP parameters.

    Returns:
        SlabResult with relaxed structure.
    """
    import numpy as np
    from ase.constraints import FixAtoms

    from ..calculator import Vasp
    from ..parameters import get_preset

    # Fix bottom layers
    if fix_bottom_layers > 0:
        positions = slab.positions[:, 2]
        sorted_z = np.sort(np.unique(np.round(positions, decimals=2)))

        if len(sorted_z) >= fix_bottom_layers:
            z_threshold = sorted_z[fix_bottom_layers - 1] + 0.1
            mask = positions <= z_threshold
            constraint = FixAtoms(mask=mask)
            slab.set_constraint(constraint)

    preset_name = 'relax-cell' if relax_cell else 'relax'
    params = get_preset(preset_name)
    params['ediffg'] = -fmax

    if dipole_correction:
        params['ldipol'] = True
        params['idipol'] = 3

    if relax_cell:
        # Only relax in-plane
        params['isif'] = 4  # Relax ions + cell shape, not volume

    params.update(calc_kwargs)

    calc = Vasp(
        atoms=slab.copy(),
        runner=runner,
        **params
    )

    calc.atoms.calc = calc
    energy = calc.atoms.get_potential_energy()
    calc.read_results()

    return SlabResult(
        atoms=calc.atoms.copy(),
        energy=energy,
        forces=calc.results.get('forces'),
        stress=calc.results.get('stress'),
        fermi_level=calc.results.get('fermi_level'),
        parameters=calc.parameters.copy(),
        directory=calc.directory,
        miller_indices=slab.info.get('miller_indices'),
        layers=slab.info.get('layers', 0),
    )


@subflow
def bulk_to_slabs_subflow(
    bulk_atoms: Atoms,
    runner: Runner | None = None,
    miller_indices: list[tuple[int, int, int]] | None = None,
    min_slab_size: float = 10.0,
    vacuum: float = 15.0,
    relax: bool = True,
    **calc_kwargs
) -> list[SlabResult]:
    """Generate and calculate multiple slabs from bulk.

    Args:
        bulk_atoms: Bulk atomic structure.
        runner: Job runner.
        miller_indices: Miller indices to generate.
        min_slab_size: Minimum slab thickness in Å.
        vacuum: Vacuum thickness in Å.
        relax: If True, relax slabs; else static only.
        **calc_kwargs: Additional VASP parameters.

    Returns:
        List of SlabResults.
    """
    slabs = make_slabs_from_bulk(
        bulk_atoms,
        miller_indices=miller_indices,
        min_slab_size=min_slab_size,
        min_vacuum_size=vacuum,
    )

    results = []
    for slab in slabs:
        if relax:
            result = slab_relax_job(slab, runner=runner, **calc_kwargs)
        else:
            result = slab_static_job(slab, runner=runner, **calc_kwargs)
        results.append(result)

    return results


@flow
def bulk_to_slabs_flow(
    bulk_atoms: Atoms,
    runner: Runner | None = None,
    miller_indices: list[tuple[int, int, int]] | None = None,
    relax_bulk: bool = True,
    **calc_kwargs
) -> dict[str, Any]:
    """Complete workflow: relax bulk, then generate and relax slabs.

    Args:
        bulk_atoms: Bulk atomic structure.
        runner: Job runner.
        miller_indices: Miller indices to generate.
        relax_bulk: If True, relax bulk first.
        **calc_kwargs: Additional VASP parameters.

    Returns:
        Dict with 'bulk' and 'slabs' results.
    """
    # Optionally relax bulk
    if relax_bulk:
        bulk_result = relax_job(
            bulk_atoms,
            runner=runner,
            relax_cell=True,
            **calc_kwargs
        )
        bulk_atoms = bulk_result.atoms
    else:
        bulk_result = static_job(bulk_atoms, runner=runner, **calc_kwargs)

    # Generate and calculate slabs
    slab_results = bulk_to_slabs_subflow(
        bulk_atoms,
        runner=runner,
        miller_indices=miller_indices,
        **calc_kwargs
    )

    return {
        'bulk': bulk_result,
        'slabs': slab_results,
    }


def calculate_surface_energy(
    slab_result: SlabResult,
    bulk_energy_per_atom: float,
) -> float:
    """Calculate surface energy from slab and bulk energies.

    Args:
        slab_result: Result from slab calculation.
        bulk_energy_per_atom: Energy per atom in bulk.

    Returns:
        Surface energy in J/m².
    """
    import numpy as np

    slab = slab_result.atoms
    n_atoms = len(slab)
    e_slab = slab_result.energy

    # Surface area (both sides)
    cell = slab.get_cell()
    area = np.linalg.norm(np.cross(cell[0], cell[1]))
    area *= 2  # Two surfaces

    # Surface energy
    e_surf_ev = (e_slab - n_atoms * bulk_energy_per_atom) / area

    # Convert eV/Å² to J/m²
    ev_per_angstrom2_to_j_per_m2 = 16.0217663

    return e_surf_ev * ev_per_angstrom2_to_j_per_m2
