"""Phonon calculation recipes with Phonopy integration.

Provides workflows for:
- Force constants calculation
- Phonon band structure
- Phonon DOS
- Thermodynamic properties
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from .core import VaspResult
from .decorators import flow, job, subflow

if TYPE_CHECKING:
    from ase import Atoms

    from ..runners.base import Runner


@dataclass
class PhononResult:
    """Result container for phonon calculations.

    Attributes:
        atoms: Unit cell atoms.
        force_constants: Force constant matrix.
        frequencies: Phonon frequencies at gamma.
        band_structure: Band structure data.
        dos: Phonon DOS data.
        thermal: Thermodynamic properties.
        phonopy: Phonopy object (if available).
    """
    atoms: Atoms | None = None
    supercell_matrix: tuple | None = None
    force_constants: Any | None = None
    frequencies: Any | None = None  # At gamma
    band_structure: dict | None = None
    dos: dict | None = None
    thermal: dict | None = None
    phonopy: Any | None = None  # Phonopy object
    directory: str = ''


def _check_phonopy():
    """Check if phonopy is available."""
    try:
        import phonopy  # noqa: F401
        return True
    except ImportError:
        return False


def get_phonopy_supercells(
    atoms: Atoms,
    supercell_matrix: tuple[int, int, int] = (2, 2, 2),
    displacement: float = 0.01,
) -> tuple[Any, list[Atoms]]:
    """Generate displaced supercells using Phonopy.

    Args:
        atoms: Unit cell structure.
        supercell_matrix: Supercell dimensions.
        displacement: Displacement distance in Å.

    Returns:
        Tuple of (Phonopy object, list of displaced supercells).

    Raises:
        ImportError: If phonopy is not installed.
    """
    if not _check_phonopy():
        raise ImportError("phonopy required: pip install phonopy")

    import os
    import tempfile

    import numpy as np
    from ase.io import write
    from phonopy import Phonopy
    from phonopy.interface.vasp import read_vasp

    # Convert ASE atoms to Phonopy format
    with tempfile.NamedTemporaryFile(suffix='.vasp', delete=False) as f:
        temp_path = f.name
        write(temp_path, atoms, format='vasp')

    try:
        unitcell = read_vasp(temp_path)
    finally:
        os.unlink(temp_path)

    # Create Phonopy object
    phonon = Phonopy(
        unitcell,
        supercell_matrix=np.diag(supercell_matrix),
        factor=521.471,  # THz to cm^-1
    )

    # Generate displacements
    phonon.generate_displacements(distance=displacement)
    supercells = phonon.supercells_with_displacements

    # Convert to ASE Atoms
    from ase import Atoms as ASEAtoms

    displaced_cells = []
    for sc in supercells:
        ase_atoms = ASEAtoms(
            symbols=sc.symbols,
            positions=sc.positions,
            cell=sc.cell,
            pbc=True,
        )
        displaced_cells.append(ase_atoms)

    return phonon, displaced_cells


@job
def phonon_displacement_job(
    atoms: Atoms,
    runner: Runner | None = None,
    **calc_kwargs
) -> VaspResult:
    """Calculate forces for a displaced supercell.

    Args:
        atoms: Displaced supercell structure.
        runner: Job runner.
        **calc_kwargs: Additional VASP parameters.

    Returns:
        VaspResult with forces.
    """
    from ..calculator import Vasp
    from ..parameters import get_phonon_params

    params = get_phonon_params('phonopy')
    params.update(calc_kwargs)

    calc = Vasp(
        atoms=atoms,
        runner=runner,
        **params
    )

    atoms.calc = calc
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()

    return VaspResult(
        atoms=atoms.copy(),
        energy=energy,
        forces=forces,
        parameters=calc.parameters.copy(),
        directory=calc.directory,
    )


@subflow
def phonon_forces_subflow(
    displaced_cells: list[Atoms],
    runner: Runner | None = None,
    **calc_kwargs
) -> list[VaspResult]:
    """Calculate forces for all displaced supercells.

    Args:
        displaced_cells: List of displaced supercell structures.
        runner: Job runner.
        **calc_kwargs: Additional VASP parameters.

    Returns:
        List of VaspResults with forces.
    """
    results = []
    for atoms in displaced_cells:
        result = phonon_displacement_job(atoms, runner=runner, **calc_kwargs)
        results.append(result)
    return results


@flow
def phonon_flow(
    atoms: Atoms,
    runner: Runner | None = None,
    supercell_matrix: tuple[int, int, int] = (2, 2, 2),
    displacement: float = 0.01,
    calculate_band_structure: bool = True,
    calculate_dos: bool = True,
    calculate_thermal: bool = True,
    temperature_range: tuple[float, float, float] = (0, 1000, 10),
    **calc_kwargs
) -> PhononResult:
    """Complete phonon calculation workflow.

    Uses Phonopy for displacement generation and post-processing.

    Args:
        atoms: Unit cell structure.
        runner: Job runner.
        supercell_matrix: Supercell dimensions (e.g., (2, 2, 2)).
        displacement: Displacement distance in Å.
        calculate_band_structure: Calculate phonon bands.
        calculate_dos: Calculate phonon DOS.
        calculate_thermal: Calculate thermal properties.
        temperature_range: (T_min, T_max, T_step) in K.
        **calc_kwargs: Additional VASP parameters.

    Returns:
        PhononResult with frequencies and properties.

    Example:
        >>> result = phonon_flow(atoms, supercell_matrix=(3, 3, 3))
        >>> print(result.frequencies)
    """
    import numpy as np

    # Generate displaced supercells
    phonon, displaced_cells = get_phonopy_supercells(
        atoms,
        supercell_matrix=supercell_matrix,
        displacement=displacement,
    )

    # Calculate forces for all displacements
    force_results = phonon_forces_subflow(
        displaced_cells,
        runner=runner,
        **calc_kwargs
    )

    # Collect forces
    forces_list = [r.forces for r in force_results]
    phonon.forces = np.array(forces_list)

    # Produce force constants
    phonon.produce_force_constants()

    result = PhononResult(
        atoms=atoms.copy(),
        supercell_matrix=supercell_matrix,
        force_constants=phonon.force_constants,
        phonopy=phonon,
        directory=force_results[0].directory if force_results else '',
    )

    # Calculate frequencies at gamma
    phonon.run_qpoints([[0, 0, 0]])
    result.frequencies = phonon.get_qpoints_dict()['frequencies'][0]

    # Band structure
    if calculate_band_structure:
        # Use automatic band path
        try:
            phonon.auto_band_structure(plot=False)
            band_dict = phonon.get_band_structure_dict()
            result.band_structure = {
                'qpoints': band_dict['qpoints'],
                'frequencies': band_dict['frequencies'],
                'labels': band_dict.get('labels', []),
            }
        except Exception:
            pass  # Skip if band structure fails

    # DOS
    if calculate_dos:
        try:
            phonon.run_mesh([20, 20, 20])
            phonon.run_total_dos()
            dos_dict = phonon.get_total_dos_dict()
            result.dos = {
                'frequency_points': dos_dict['frequency_points'],
                'total_dos': dos_dict['total_dos'],
            }
        except Exception:
            pass

    # Thermal properties
    if calculate_thermal:
        try:
            phonon.run_mesh([20, 20, 20])
            t_min, t_max, t_step = temperature_range
            phonon.run_thermal_properties(
                t_min=t_min,
                t_max=t_max,
                t_step=t_step,
            )
            tp_dict = phonon.get_thermal_properties_dict()
            result.thermal = {
                'temperatures': tp_dict['temperatures'],
                'free_energy': tp_dict['free_energy'],
                'entropy': tp_dict['entropy'],
                'heat_capacity': tp_dict['heat_capacity'],
            }
        except Exception:
            pass

    return result


@job
def dfpt_phonon_job(
    atoms: Atoms,
    runner: Runner | None = None,
    **calc_kwargs
) -> PhononResult:
    """DFPT phonon calculation using VASP's internal routines.

    Uses IBRION=8 for symmetry-reduced DFPT.
    Only calculates gamma-point phonons (q=0).

    Args:
        atoms: Unit cell structure.
        runner: Job runner.
        **calc_kwargs: Additional VASP parameters.

    Returns:
        PhononResult with gamma frequencies.
    """
    from ..calculator import Vasp
    from ..parameters import get_phonon_params

    params = get_phonon_params('dfpt')
    params.update(calc_kwargs)

    calc = Vasp(
        atoms=atoms,
        runner=runner,
        **params
    )

    atoms.calc = calc
    atoms.get_potential_energy()

    # Parse frequencies from OUTCAR
    frequencies = _parse_dfpt_frequencies(calc.directory)

    return PhononResult(
        atoms=atoms.copy(),
        frequencies=frequencies,
        directory=calc.directory,
    )


def _parse_dfpt_frequencies(directory: str) -> list[float]:
    """Parse phonon frequencies from DFPT OUTCAR.

    Args:
        directory: Calculation directory.

    Returns:
        List of frequencies in cm^-1.
    """
    import os
    import re

    outcar_path = os.path.join(directory, 'OUTCAR')
    if not os.path.exists(outcar_path):
        return []

    frequencies = []
    with open(outcar_path) as f:
        content = f.read()

    # Pattern for DFPT frequencies
    # "f  =    5.123456 THz    32.18765 2PiTHz  170.87654 cm-1    21.18765 meV"
    pattern = r'f\s*[/=i]*\s*([-\d.]+)\s+THz.*?([-\d.]+)\s+cm-1'
    matches = re.findall(pattern, content)

    for _thz, cm1 in matches:
        frequencies.append(float(cm1))

    return frequencies
