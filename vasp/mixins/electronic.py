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

        if "fermi_level" in self.results:
            return self.results["fermi_level"]

        # Try to read from vasprun.xml
        return self._read_fermi_from_vasprun()

    def _read_fermi_from_vasprun(self) -> float:
        """Read Fermi level from vasprun.xml."""
        from xml.etree import ElementTree

        vasprun = os.path.join(self.directory, "vasprun.xml")
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

        vasprun = os.path.join(self.directory, "vasprun.xml")
        tree = ElementTree.parse(vasprun)

        eigenvalues = []
        occupations = []

        # Find all spin blocks
        for spin_set in tree.findall(".//eigenvalues/array/set/set"):
            spin_eigs = []
            spin_occs = []

            for kpt_set in spin_set.findall("set"):
                kpt_eigs = []
                kpt_occs = []

                for r in kpt_set.findall("r"):
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

        vasprun = os.path.join(self.directory, "vasprun.xml")
        tree = ElementTree.parse(vasprun)

        kpts_elem = tree.find(".//kpoints/varray[@name='kpointlist']")
        kpts = []
        for v in kpts_elem.findall("v"):
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

        vasprun = os.path.join(self.directory, "vasprun.xml")
        tree = ElementTree.parse(vasprun)

        weights_elem = tree.find(".//kpoints/varray[@name='weights']")
        weights = []
        for v in weights_elem.findall("v"):
            weights.append(float(v.text))

        return np.array(weights)

    def get_number_of_spins(self) -> int:
        """Get number of spin channels.

        Returns:
            1 for non-spin-polarized, 2 for spin-polarized.
        """
        # First check parameters
        if "ispin" in self.parameters:
            return self.parameters["ispin"]

        # Try to read from INCAR
        incar = os.path.join(self.directory, "INCAR")
        if os.path.exists(incar):
            with open(incar) as f:
                for line in f:
                    if "ISPIN" in line.upper():
                        parts = line.split("=")
                        if len(parts) >= 2:
                            try:
                                return int(parts[1].split()[0])
                            except (ValueError, IndexError):
                                pass

        # Try vasprun.xml
        vasprun = os.path.join(self.directory, "vasprun.xml")
        if os.path.exists(vasprun):
            from xml.etree import ElementTree

            tree = ElementTree.parse(vasprun)
            ispin_elem = tree.find(
                './/parameters/separator[@name="electronic spin"]/i[@name="ISPIN"]'
            )
            if ispin_elem is not None:
                return int(ispin_elem.text.strip())

        return 1  # Default

    def get_dos(
        self, spin: int | None = None, efermi: float | None = None
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

        vasprun = os.path.join(self.directory, "vasprun.xml")
        tree = ElementTree.parse(vasprun)

        if efermi is None:
            efermi = self.get_fermi_level()

        dos_elem = tree.find(".//dos/total/array/set")
        if dos_elem is None:
            raise ValueError("DOS not found in vasprun.xml")

        # Read DOS data
        energies = []
        dos_data = []

        for spin_set in dos_elem.findall("set"):
            spin_dos = []
            for r in spin_set.findall("r"):
                parts = r.text.split()
                if not energies or len(energies) < len(spin_set.findall("r")):
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

    def get_ados(
        self, atom_index: int, orbital: str, spin: int | None = None, efermi: float | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        """Get atom-projected density of states.

        Requires LORBIT=10, 11, or 12 in the calculation.

        Args:
            atom_index: Atom index (0-based).
            orbital: Orbital type ('s', 'p', 'd', 'f', or 'total').
            spin: Spin channel (None for total, 0 or 1 for spin-polarized).
            efermi: Fermi level for shifting energies (default: use calculated).

        Returns:
            Tuple of (energies, pdos) arrays.
        """
        self.update()
        from xml.etree import ElementTree

        vasprun = os.path.join(self.directory, "vasprun.xml")
        tree = ElementTree.parse(vasprun)

        if efermi is None:
            efermi = self.get_fermi_level()

        # Find projected DOS
        partial = tree.find(".//dos/partial")
        if partial is None:
            raise ValueError(
                "Projected DOS not found. Set LORBIT=10, 11, or 12 in your calculation."
            )

        # Get orbital names from the field definitions
        fields = partial.find(".//set/set/set/r")
        if fields is None:
            raise ValueError("Could not parse projected DOS structure")

        # Determine orbital columns from header
        # VASP typically outputs: energy s py pz px dxy dyz dz2 dxz dx2 (for LORBIT=11)
        # or: energy s p d (for LORBIT=10)
        orbital_map = {
            "s": [1],
            "p": [2, 3, 4] if len(fields.text.split()) > 5 else [2],
            "d": [5, 6, 7, 8, 9] if len(fields.text.split()) > 5 else [3],
            "f": [10, 11, 12, 13, 14, 15, 16] if len(fields.text.split()) > 10 else [4],
        }

        # Find the ion set for the requested atom
        ions = partial.findall('.//set[@comment="ion"]/set')
        if not ions:
            # Try alternative structure
            ions = partial.findall(".//array/set/set")

        if atom_index >= len(ions):
            raise IndexError(f"Atom index {atom_index} out of range (have {len(ions)} atoms)")

        ion_set = ions[atom_index]

        # Read DOS for each spin
        energies = []
        pdos_data = []

        for spin_set in ion_set.findall("set"):
            spin_pdos = []
            for r in spin_set.findall("r"):
                parts = [float(x) for x in r.text.split()]
                if not energies or len(energies) < len(spin_set.findall("r")):
                    energies.append(parts[0])

                # Sum requested orbital columns
                if orbital.lower() == "total":
                    # Sum all orbitals (skip energy column)
                    spin_pdos.append(sum(parts[1:]))
                else:
                    cols = orbital_map.get(orbital.lower(), [1])
                    # Handle case where orbital columns don't exist
                    valid_cols = [c for c in cols if c < len(parts)]
                    if valid_cols:
                        spin_pdos.append(sum(parts[c] for c in valid_cols))
                    else:
                        spin_pdos.append(0.0)

            pdos_data.append(np.array(spin_pdos))

        energies = np.array(energies) - efermi

        if spin is not None:
            return energies, pdos_data[spin]
        elif len(pdos_data) == 2:
            return energies, pdos_data[0] + pdos_data[1]
        else:
            return energies, pdos_data[0]

    def get_magnetic_moment(self) -> float:
        """Get total magnetic moment.

        Returns:
            Total magnetic moment in Bohr magnetons.
        """
        self.update()

        if "magmom" in self.results:
            return self.results["magmom"]

        # Try to read from OUTCAR
        outcar = os.path.join(self.directory, "OUTCAR")
        if os.path.exists(outcar):
            with open(outcar) as f:
                content = f.read()
            match = re.search(r"number of electron\s+[\d.]+\s+magnetization\s+([-\d.]+)", content)
            if match:
                return float(match.group(1))

        return 0.0

    def get_magnetic_moments(self) -> np.ndarray:
        """Get per-atom magnetic moments.

        Returns:
            1D array of magnetic moments per atom.
        """
        self.update()

        outcar = os.path.join(self.directory, "OUTCAR")
        if not os.path.exists(outcar):
            raise FileNotFoundError(f"OUTCAR not found in {self.directory}")

        with open(outcar) as f:
            content = f.read()

        # Find magnetization section
        pattern = r"magnetization \(x\)\s*\n.*?\n((?:\s*\d+\s+[-\d.]+\s*\n)+)"
        matches = re.findall(pattern, content, re.DOTALL)

        if matches:
            magmoms = []
            for line in matches[-1].strip().split("\n"):
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

        homo = float("-inf")
        lumo = float("inf")

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
