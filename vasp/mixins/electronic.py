"""Electronic properties mixin for VASP calculator.

Provides methods for accessing electronic structure data:
- Eigenvalues and occupation numbers
- Fermi level
- Density of states (DOS)
- Band structure data
- Charge density information
"""

from __future__ import annotations

import os
import re
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    pass


class ElectronicMixin:
    """Mixin providing electronic structure property methods.

    Requires:
    - directory: str - Calculation directory path
    - results: dict - Results dictionary
    - update(): method - Ensure calculation is complete
    """

    def get_fermi_level(self) -> float:
        """Get the Fermi level in eV.

        Returns:
            Fermi energy in eV.

        Raises:
            KeyError: If Fermi level not found in results.
        """
        self.update()

        if 'fermi_level' in self.results:
            return self.results['fermi_level']

        # Try to read from vasprun.xml
        return self._read_fermi_from_vasprun()

    def _read_fermi_from_vasprun(self) -> float:
        """Read Fermi level from vasprun.xml."""
        from xml.etree import ElementTree

        vasprun = os.path.join(self.directory, 'vasprun.xml')
        if not os.path.exists(vasprun):
            raise FileNotFoundError(f"vasprun.xml not found in {self.directory}")

        tree = ElementTree.parse(vasprun)
        fermi = tree.find('.//dos/i[@name="efermi"]')
        if fermi is not None:
            return float(fermi.text)

        raise KeyError("Fermi level not found in vasprun.xml")

    def get_eigenvalues(self, kpt: int = 0, spin: int = 0) -> np.ndarray:
        """Get eigenvalues for a specific k-point and spin.

        Args:
            kpt: K-point index (0-based).
            spin: Spin index (0 or 1).

        Returns:
            1D array of eigenvalues in eV.
        """
        self.update()
        eigenvalues = self._read_eigenvalues()

        if spin >= len(eigenvalues):
            raise IndexError(f"Spin index {spin} out of range")
        if kpt >= len(eigenvalues[spin]):
            raise IndexError(f"K-point index {kpt} out of range")

        return eigenvalues[spin][kpt]

    def get_occupation_numbers(self, kpt: int = 0, spin: int = 0) -> np.ndarray:
        """Get occupation numbers for a specific k-point and spin.

        Args:
            kpt: K-point index (0-based).
            spin: Spin index (0 or 1).

        Returns:
            1D array of occupation numbers.
        """
        self.update()
        _, occupations = self._read_eigenvalues_and_occupations()

        if spin >= len(occupations):
            raise IndexError(f"Spin index {spin} out of range")
        if kpt >= len(occupations[spin]):
            raise IndexError(f"K-point index {kpt} out of range")

        return occupations[spin][kpt]

    def _read_eigenvalues(self) -> list[list[np.ndarray]]:
        """Read eigenvalues from vasprun.xml.

        Returns:
            Nested list: eigenvalues[spin][kpt] = 1D array
        """
        eigenvalues, _ = self._read_eigenvalues_and_occupations()
        return eigenvalues

    def _read_eigenvalues_and_occupations(self) -> tuple[list, list]:
        """Read eigenvalues and occupations from vasprun.xml.

        Returns:
            Tuple of (eigenvalues, occupations) as nested lists.
        """
        from xml.etree import ElementTree

        vasprun = os.path.join(self.directory, 'vasprun.xml')
        tree = ElementTree.parse(vasprun)

        eigenvalues = []
        occupations = []

        # Find all spin blocks
        for spin_set in tree.findall('.//eigenvalues/array/set/set'):
            spin_eigs = []
            spin_occs = []

            for kpt_set in spin_set.findall('set'):
                kpt_eigs = []
                kpt_occs = []

                for r in kpt_set.findall('r'):
                    parts = r.text.split()
                    kpt_eigs.append(float(parts[0]))
                    kpt_occs.append(float(parts[1]))

                spin_eigs.append(np.array(kpt_eigs))
                spin_occs.append(np.array(kpt_occs))

            eigenvalues.append(spin_eigs)
            occupations.append(spin_occs)

        return eigenvalues, occupations

    def get_ibz_k_points(self, cartesian: bool = True) -> np.ndarray:
        """Get irreducible Brillouin zone k-points.

        Args:
            cartesian: If True, return in Cartesian coordinates.
                If False, return in reciprocal lattice coordinates.

        Returns:
            Array of k-points with shape (nkpts, 3).
        """
        self.update()
        from xml.etree import ElementTree

        vasprun = os.path.join(self.directory, 'vasprun.xml')
        tree = ElementTree.parse(vasprun)

        kpts_elem = tree.find(".//kpoints/varray[@name='kpointlist']")
        kpts = []
        for v in kpts_elem.findall('v'):
            kpts.append([float(x) for x in v.text.split()])

        kpts = np.array(kpts)

        if cartesian:
            # Transform to Cartesian
            cell = self.atoms.get_cell()
            recip = 2 * np.pi * np.linalg.inv(cell).T
            kpts = kpts @ recip

        return kpts

    def get_k_point_weights(self) -> np.ndarray:
        """Get k-point weights.

        Returns:
            1D array of weights summing to 1.
        """
        self.update()
        from xml.etree import ElementTree

        vasprun = os.path.join(self.directory, 'vasprun.xml')
        tree = ElementTree.parse(vasprun)

        weights_elem = tree.find(".//kpoints/varray[@name='weights']")
        weights = []
        for v in weights_elem.findall('v'):
            weights.append(float(v.text))

        return np.array(weights)

    def get_number_of_spins(self) -> int:
        """Get number of spin channels.

        Returns:
            1 for non-spin-polarized, 2 for spin-polarized.
        """
        self.update()
        ispin = self.parameters.get('ispin', 1)
        return ispin

    def get_dos(
        self,
        spin: int | None = None,
        efermi: float | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        """Get density of states.

        Args:
            spin: Spin channel (None for total, 0 or 1 for spin-polarized).
            efermi: Fermi level for shifting energies (default: use calculated).

        Returns:
            Tuple of (energies, dos) arrays.
        """
        self.update()
        from xml.etree import ElementTree

        vasprun = os.path.join(self.directory, 'vasprun.xml')
        tree = ElementTree.parse(vasprun)

        if efermi is None:
            efermi = self.get_fermi_level()

        dos_elem = tree.find('.//dos/total/array/set')
        if dos_elem is None:
            raise ValueError("DOS not found in vasprun.xml")

        # Read DOS data
        energies = []
        dos_data = []

        for spin_set in dos_elem.findall('set'):
            spin_dos = []
            for r in spin_set.findall('r'):
                parts = r.text.split()
                if not energies or len(energies) < len(spin_set.findall('r')):
                    energies.append(float(parts[0]))
                spin_dos.append(float(parts[1]))
            dos_data.append(np.array(spin_dos))

        energies = np.array(energies) - efermi

        if spin is not None:
            return energies, dos_data[spin]
        elif len(dos_data) == 2:
            return energies, dos_data[0] + dos_data[1]
        else:
            return energies, dos_data[0]

    def get_magnetic_moment(self) -> float:
        """Get total magnetic moment.

        Returns:
            Total magnetic moment in Bohr magnetons.
        """
        self.update()

        if 'magmom' in self.results:
            return self.results['magmom']

        # Try to read from OUTCAR
        outcar = os.path.join(self.directory, 'OUTCAR')
        if os.path.exists(outcar):
            with open(outcar) as f:
                content = f.read()
            match = re.search(r'number of electron\s+[\d.]+\s+magnetization\s+([-\d.]+)', content)
            if match:
                return float(match.group(1))

        return 0.0

    def get_magnetic_moments(self) -> np.ndarray:
        """Get per-atom magnetic moments.

        Returns:
            1D array of magnetic moments per atom.
        """
        self.update()

        outcar = os.path.join(self.directory, 'OUTCAR')
        if not os.path.exists(outcar):
            raise FileNotFoundError(f"OUTCAR not found in {self.directory}")

        with open(outcar) as f:
            content = f.read()

        # Find magnetization section
        pattern = r'magnetization \(x\)\s*\n.*?\n((?:\s*\d+\s+[-\d.]+\s*\n)+)'
        matches = re.findall(pattern, content, re.DOTALL)

        if matches:
            magmoms = []
            for line in matches[-1].strip().split('\n'):
                parts = line.split()
                if len(parts) >= 2:
                    magmoms.append(float(parts[1]))
            return np.array(magmoms)[self.resort]

        return np.zeros(len(self.atoms))

    def get_homo_lumo(self) -> tuple[float, float]:
        """Get HOMO and LUMO energies.

        Returns:
            Tuple of (HOMO, LUMO) energies in eV relative to vacuum.
        """
        self.update()

        efermi = self.get_fermi_level()
        eigenvalues = self._read_eigenvalues()

        homo = float('-inf')
        lumo = float('inf')

        for spin_eigs in eigenvalues:
            for kpt_eigs in spin_eigs:
                # Find highest occupied and lowest unoccupied
                below_fermi = kpt_eigs[kpt_eigs <= efermi + 0.01]
                above_fermi = kpt_eigs[kpt_eigs > efermi + 0.01]

                if len(below_fermi) > 0:
                    homo = max(homo, below_fermi.max())
                if len(above_fermi) > 0:
                    lumo = min(lumo, above_fermi.min())

        return homo, lumo

    def get_band_gap(self) -> float:
        """Get the band gap.

        Returns:
            Band gap in eV (0 for metals).
        """
        homo, lumo = self.get_homo_lumo()

        if lumo > homo:
            return lumo - homo
        return 0.0
