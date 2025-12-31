"""Dynamics mixin for NEB and vibrational calculations.

Provides methods for:
- Nudged Elastic Band (NEB) calculations
- Vibrational frequency analysis
- Phonon calculations
"""

from __future__ import annotations

import os
import re
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ase import Atoms


class DynamicsMixin:
    """Mixin providing NEB and vibrational analysis methods.

    Requires:
    - directory: str - Calculation directory path
    - atoms: Atoms - ASE Atoms object
    - parameters: dict - VASP parameters
    - update(): method - Ensure calculation is complete
    """

    def get_vibrational_frequencies(self) -> np.ndarray:
        """Get vibrational frequencies from OUTCAR.

        Requires IBRION=5, 6, 7, or 8 calculation.

        Returns:
            1D array of frequencies in cm^-1.
            Imaginary frequencies are returned as negative values.
        """
        self.update()
        frequencies, _ = self._parse_vibrations()
        return frequencies

    def get_vibrational_modes(self, massweighted: bool = False) -> tuple[np.ndarray, np.ndarray]:
        """Get vibrational frequencies and eigenvectors.

        Args:
            massweighted: If True, return mass-weighted eigenvectors.

        Returns:
            Tuple of (frequencies, eigenvectors) where:
            - frequencies: 1D array in cm^-1
            - eigenvectors: 3D array of shape (nmodes, natoms, 3)
        """
        self.update()
        frequencies, eigenvectors = self._parse_vibrations()

        if not massweighted:
            # Un-mass-weight the eigenvectors
            masses = self.atoms.get_masses()
            for i in range(len(frequencies)):
                for j in range(len(self.atoms)):
                    eigenvectors[i, j] *= np.sqrt(masses[j])
                # Normalize
                norm = np.linalg.norm(eigenvectors[i])
                if norm > 0:
                    eigenvectors[i] /= norm

        return frequencies, eigenvectors

    def _parse_vibrations(self) -> tuple[np.ndarray, np.ndarray]:
        """Parse vibrational data from OUTCAR.

        Returns:
            Tuple of (frequencies, eigenvectors).
        """
        outcar = os.path.join(self.directory, "OUTCAR")

        if not os.path.exists(outcar):
            raise FileNotFoundError(f"OUTCAR not found in {self.directory}")

        with open(outcar) as f:
            content = f.read()

        # Check for vibrational calculation
        if "Eigenvectors and eigenvalues" not in content:
            raise ValueError(
                "No vibrational data found. " "Use IBRION=5, 6, 7, or 8 for vibrations."
            )

        frequencies = []
        eigenvectors = []
        natoms = len(self.atoms)

        # Parse each mode
        pattern = (
            r"(\d+)\s+f\s*(/i)?=\s*([\d.]+)\s+THz\s+"
            r"([\d.]+)\s+2PiTHz\s+"
            r"([\d.]+)\s+cm-1\s+"
            r"([\d.]+)\s+meV"
        )

        re.split(r"\d+\s+f\s*/?\s*i?=", content)[1:]

        for match in re.finditer(pattern, content):
            int(match.group(1))
            is_imaginary = match.group(2) == "/i"
            freq_cm = float(match.group(5))

            if is_imaginary:
                freq_cm = -freq_cm

            frequencies.append(freq_cm)

        # Parse eigenvectors

        # Simpler approach: find all displacement blocks
        blocks = re.findall(
            r"X\s+Y\s+Z\s+dx\s+dy\s+dz\s*\n((?:\s*[\d.-]+\s+[\d.-]+\s+[\d.-]+\s+[\d.-]+\s+[\d.-]+\s+[\d.-]+\s*\n)+)",
            content,
        )

        for block in blocks:
            mode_vecs = []
            for line in block.strip().split("\n"):
                parts = line.split()
                if len(parts) >= 6:
                    dx, dy, dz = float(parts[3]), float(parts[4]), float(parts[5])
                    mode_vecs.append([dx, dy, dz])

            if len(mode_vecs) == natoms:
                eigenvectors.append(mode_vecs)

        frequencies = np.array(frequencies)
        eigenvectors = (
            np.array(eigenvectors) if eigenvectors else np.zeros((len(frequencies), natoms, 3))
        )

        return frequencies, eigenvectors

    def get_zero_point_energy(self) -> float:
        """Calculate zero-point energy from vibrational frequencies.

        Returns:
            Zero-point energy in eV.
        """
        frequencies = self.get_vibrational_frequencies()

        # Only use real (positive) frequencies
        real_freqs = frequencies[frequencies > 0]

        # ZPE = 0.5 * sum(hbar * omega)
        # Convert cm^-1 to eV: 1 cm^-1 = 1.239841984e-4 eV
        cm_to_ev = 1.239841984e-4 / 2.0  # includes 0.5 factor

        zpe = np.sum(real_freqs) * cm_to_ev

        return zpe

    def get_infrared_intensities(self) -> np.ndarray:
        """Get IR intensities from OUTCAR.

        Requires LEPSILON=.TRUE. or LCALCEPS=.TRUE.

        Returns:
            1D array of IR intensities.
        """
        self.update()

        outcar = os.path.join(self.directory, "OUTCAR")
        with open(outcar) as f:
            content = f.read()

        # Look for Born effective charges and compute intensities
        if "BORN EFFECTIVE CHARGES" not in content:
            raise ValueError(
                "Born effective charges not found. " "Use LEPSILON=.TRUE. for IR intensities."
            )

        # Parse Born effective charges
        pattern = r"ion\s+\d+\s*\n((?:\s*\d+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n){3})"
        matches = re.findall(pattern, content)

        if not matches:
            raise ValueError("Could not parse Born effective charges")

        born_charges = []
        for match in matches:
            ion_charges = []
            for line in match.strip().split("\n"):
                parts = line.split()
                if len(parts) >= 4:
                    ion_charges.append([float(parts[1]), float(parts[2]), float(parts[3])])
            born_charges.append(ion_charges)

        born_charges = np.array(born_charges)

        # Get vibrational modes
        frequencies, eigenvectors = self.get_vibrational_modes()

        # Calculate IR intensities
        # I = |sum_i Z*_i . e_i|^2
        intensities = []
        for mode in eigenvectors:
            dipole = np.zeros(3)
            for i, evec in enumerate(mode):
                if i < len(born_charges):
                    dipole += born_charges[i] @ evec
            intensity = np.linalg.norm(dipole) ** 2
            intensities.append(intensity)

        return np.array(intensities)

    def get_neb_images(self) -> list[Atoms]:
        """Get NEB images from calculation.

        Returns:
            List of Atoms objects for each NEB image.
        """
        from ase.io import read

        images = []

        # Check for numbered directories (00, 01, 02, ...)
        for i in range(100):
            image_dir = os.path.join(self.directory, f"{i:02d}")
            if not os.path.exists(image_dir):
                break

            contcar = os.path.join(image_dir, "CONTCAR")
            poscar = os.path.join(image_dir, "POSCAR")

            if os.path.exists(contcar) and os.path.getsize(contcar) > 0:
                atoms = read(contcar, format="vasp")
            elif os.path.exists(poscar):
                atoms = read(poscar, format="vasp")
            else:
                continue

            # Try to get energy
            outcar = os.path.join(image_dir, "OUTCAR")
            if os.path.exists(outcar):
                with open(outcar) as f:
                    content = f.read()
                match = re.search(r"free  energy   TOTEN\s*=\s*([-\d.]+)", content)
                if match:
                    atoms.info["energy"] = float(match.group(1))

            images.append(atoms)

        if not images:
            raise FileNotFoundError("No NEB images found")

        return images

    def get_neb_barrier(self) -> tuple[float, float]:
        """Get NEB activation energy barriers.

        Returns:
            Tuple of (forward_barrier, reverse_barrier) in eV.
        """
        images = self.get_neb_images()

        energies = []
        for img in images:
            if "energy" in img.info:
                energies.append(img.info["energy"])
            else:
                raise ValueError("Energy not found for all NEB images")

        energies = np.array(energies)

        # Find transition state (maximum)
        ts_idx = np.argmax(energies)
        ts_energy = energies[ts_idx]

        # Forward and reverse barriers
        forward = ts_energy - energies[0]
        reverse = ts_energy - energies[-1]

        return forward, reverse

    def get_neb_path(self) -> tuple[np.ndarray, np.ndarray]:
        """Get NEB reaction path data.

        Returns:
            Tuple of (reaction_coordinate, energies) arrays.
        """
        images = self.get_neb_images()

        # Calculate reaction coordinate as cumulative distance
        distances = [0.0]
        for i in range(1, len(images)):
            d = np.linalg.norm(images[i].positions - images[i - 1].positions)
            distances.append(distances[-1] + d)

        # Normalize to [0, 1]
        distances = np.array(distances)
        distances /= distances[-1]

        energies = np.array([img.info.get("energy", np.nan) for img in images])

        # Reference to initial state
        energies -= energies[0]

        return distances, energies
