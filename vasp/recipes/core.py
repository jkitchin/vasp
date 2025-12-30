"""Core VASP calculation recipes.

Provides basic jobs for static and relaxation calculations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from .decorators import flow, job

if TYPE_CHECKING:
    from ase import Atoms

    from ..runners.base import Runner


@dataclass
class VaspResult:
    """Result container for VASP calculations.

    Attributes:
        atoms: Final atomic structure.
        energy: Total energy in eV.
        forces: Forces array (N, 3) in eV/Å.
        stress: Stress tensor in Voigt notation.
        fermi_level: Fermi energy in eV.
        band_gap: Band gap in eV (if calculated).
        parameters: VASP parameters used.
        directory: Calculation directory.
    """
    atoms: Atoms | None = None
    energy: float | None = None
    forces: Any | None = None
    stress: Any | None = None
    fermi_level: float | None = None
    band_gap: float | None = None
    parameters: dict = field(default_factory=dict)
    directory: str = ''
    converged: bool = True
    nsteps: int = 0


@job
def static_job(
    atoms: Atoms,
    runner: Runner | None = None,
    preset: str | None = None,
    copy_files: list[str] | None = None,
    **calc_kwargs
) -> VaspResult:
    """Run a static (single-point) VASP calculation.

    Args:
        atoms: Input atomic structure.
        runner: Job runner (LocalRunner, SlurmRunner, etc.).
        preset: Parameter preset name.
        copy_files: Files to copy from previous calculation.
        **calc_kwargs: Additional VASP parameters.

    Returns:
        VaspResult with energy, forces, stress.

    Example:
        >>> from ase.build import bulk
        >>> si = bulk('Si')
        >>> result = static_job(si, encut=400, kpts=(4,4,4))
        >>> print(result.energy)
    """
    from ..calculator import Vasp
    from ..parameters import get_preset

    # Build parameters
    params = {}
    if preset:
        params.update(get_preset(preset))
    else:
        params.update(get_preset('static'))
    params.update(calc_kwargs)

    # Create calculator
    calc = Vasp(
        atoms=atoms,
        runner=runner,
        **params
    )

    # Run calculation
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    stress = atoms.get_stress()

    # Gather results
    return VaspResult(
        atoms=atoms.copy(),
        energy=energy,
        forces=forces,
        stress=stress,
        fermi_level=calc.results.get('fermi_level'),
        parameters=calc.parameters.copy(),
        directory=calc.directory,
    )


@job
def relax_job(
    atoms: Atoms,
    runner: Runner | None = None,
    relax_cell: bool = False,
    fmax: float = 0.02,
    steps: int = 100,
    preset: str | None = None,
    **calc_kwargs
) -> VaspResult:
    """Run a VASP relaxation calculation.

    Uses VASP's internal optimizer (IBRION=2 conjugate gradient).

    Args:
        atoms: Input atomic structure.
        runner: Job runner.
        relax_cell: If True, also relax cell shape/volume (ISIF=3).
        fmax: Force convergence criterion in eV/Å.
        steps: Maximum ionic steps.
        preset: Parameter preset name.
        **calc_kwargs: Additional VASP parameters.

    Returns:
        VaspResult with relaxed structure.

    Example:
        >>> result = relax_job(atoms, relax_cell=True, fmax=0.01)
        >>> relaxed_atoms = result.atoms
    """
    from ..calculator import Vasp
    from ..parameters import get_preset

    # Build parameters
    params = {}
    if preset:
        params.update(get_preset(preset))
    else:
        preset_name = 'relax-cell' if relax_cell else 'relax'
        params.update(get_preset(preset_name))

    params['nsw'] = steps
    params['ediffg'] = -fmax  # Negative for force criterion
    params.update(calc_kwargs)

    # Create calculator
    calc = Vasp(
        atoms=atoms.copy(),
        runner=runner,
        **params
    )

    # Run calculation
    calc.atoms.calc = calc
    energy = calc.atoms.get_potential_energy()

    # Read relaxed structure
    calc.read_results()
    relaxed_atoms = calc.atoms.copy()

    return VaspResult(
        atoms=relaxed_atoms,
        energy=energy,
        forces=calc.results.get('forces'),
        stress=calc.results.get('stress'),
        fermi_level=calc.results.get('fermi_level'),
        parameters=calc.parameters.copy(),
        directory=calc.directory,
        nsteps=params['nsw'],
    )


@flow
def double_relax_flow(
    atoms: Atoms,
    runner: Runner | None = None,
    relax_cell: bool = True,
    fmax: float = 0.02,
    **calc_kwargs
) -> VaspResult:
    """Two-step relaxation for better convergence.

    First relaxation uses coarse settings, second uses fine settings.

    Args:
        atoms: Input atomic structure.
        runner: Job runner.
        relax_cell: If True, also relax cell.
        fmax: Final force convergence criterion.
        **calc_kwargs: Additional VASP parameters.

    Returns:
        VaspResult from final relaxation.
    """
    # Coarse relaxation
    coarse_kwargs = calc_kwargs.copy()
    coarse_kwargs.setdefault('ediff', 1e-5)
    coarse_kwargs.setdefault('prec', 'Normal')

    result1 = relax_job(
        atoms,
        runner=runner,
        relax_cell=relax_cell,
        fmax=fmax * 2,  # Looser criterion
        **coarse_kwargs
    )

    # Fine relaxation
    fine_kwargs = calc_kwargs.copy()
    fine_kwargs.setdefault('ediff', 1e-6)
    fine_kwargs.setdefault('prec', 'Accurate')

    result2 = relax_job(
        result1.atoms,
        runner=runner,
        relax_cell=relax_cell,
        fmax=fmax,
        **fine_kwargs
    )

    return result2


@job
def static_from_relax_job(
    relax_result: VaspResult,
    runner: Runner | None = None,
    **calc_kwargs
) -> VaspResult:
    """Run static calculation on relaxed structure.

    Args:
        relax_result: Result from relax_job.
        runner: Job runner.
        **calc_kwargs: Additional VASP parameters.

    Returns:
        VaspResult from static calculation.
    """
    return static_job(
        relax_result.atoms,
        runner=runner,
        **calc_kwargs
    )


@flow
def relax_and_static_flow(
    atoms: Atoms,
    runner: Runner | None = None,
    relax_cell: bool = False,
    fmax: float = 0.02,
    **calc_kwargs
) -> tuple[VaspResult, VaspResult]:
    """Relaxation followed by static calculation.

    Args:
        atoms: Input atomic structure.
        runner: Job runner.
        relax_cell: If True, also relax cell.
        fmax: Force convergence criterion.
        **calc_kwargs: Additional VASP parameters.

    Returns:
        Tuple of (relax_result, static_result).
    """
    relax_result = relax_job(
        atoms,
        runner=runner,
        relax_cell=relax_cell,
        fmax=fmax,
        **calc_kwargs
    )

    static_result = static_job(
        relax_result.atoms,
        runner=runner,
        **calc_kwargs
    )

    return relax_result, static_result
