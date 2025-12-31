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

        # Load atoms if not already loaded
        if self.atoms is None:
            self.load_atoms()
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

    def get_neb_path(
        self,
        initial_energy: float | None = None,
        final_energy: float | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Get NEB reaction path data.

        Args:
            initial_energy: Energy of initial state (if not in NEB output).
            final_energy: Energy of final state (if not in NEB output).

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

        # Use provided endpoint energies if available
        if initial_energy is not None:
            energies[0] = initial_energy
        if final_energy is not None:
            energies[-1] = final_energy

        # Try to get endpoint energies from original NEB images
        if self.neb_images is not None:
            if np.isnan(energies[0]):
                # Check original initial image for energy
                init_img = self.neb_images[0]
                if "energy" in init_img.info:
                    energies[0] = init_img.info["energy"]
                elif init_img.calc is not None:
                    try:
                        energies[0] = init_img.get_potential_energy()
                    except Exception:
                        pass
            if np.isnan(energies[-1]):
                # Check original final image for energy
                final_img = self.neb_images[-1]
                if "energy" in final_img.info:
                    energies[-1] = final_img.info["energy"]
                elif final_img.calc is not None:
                    try:
                        energies[-1] = final_img.get_potential_energy()
                    except Exception:
                        pass

        # Reference to initial state (use first non-NaN if initial is NaN)
        ref_energy = energies[0]
        if np.isnan(ref_energy):
            # Find first non-NaN energy
            for e in energies:
                if not np.isnan(e):
                    ref_energy = e
                    break

        energies = energies - ref_energy

        return distances, energies

    def get_neb(self) -> tuple[list[Atoms], list[float | None]]:
        """Get NEB images and their energies.

        Returns:
            Tuple of (images, energies) where images is a list of Atoms
            and energies is a list of energies (None if not available).
        """
        images = self.get_neb_images()
        energies = []

        for img in images:
            energy = img.info.get("energy")
            if energy is None and img.calc is not None:
                try:
                    energy = img.get_potential_energy()
                except Exception:
                    pass
            energies.append(energy)

        return images, energies

    def plot_neb(
        self,
        ax=None,
        show: bool = True,
        label: str | None = None,
        initial_energy: float | None = None,
        final_energy: float | None = None,
        **kwargs,
    ):
        """Plot NEB energy profile.

        Args:
            ax: Matplotlib axes. If None, creates new figure.
            show: If True, call plt.show().
            label: Legend label for the curve.
            initial_energy: Energy of initial state (if not computed in NEB).
            final_energy: Energy of final state (if not computed in NEB).
            **kwargs: Passed to plt.plot().

        Returns:
            Matplotlib axes object.
        """
        import matplotlib.pyplot as plt

        distances, energies = self.get_neb_path(
            initial_energy=initial_energy, final_energy=final_energy
        )

        if ax is None:
            fig, ax = plt.subplots()

        # Filter out NaN values for spline fitting
        valid = ~np.isnan(energies)
        x_valid = distances[valid]
        y_valid = energies[valid]

        # Plot data points
        ax.plot(x_valid, y_valid, "o", markersize=8, label=label, **kwargs)

        # Add spline fit if we have enough points
        if len(x_valid) >= 4:
            from scipy.interpolate import CubicSpline

            # Create smooth spline through points
            cs = CubicSpline(x_valid, y_valid)
            x_smooth = np.linspace(x_valid.min(), x_valid.max(), 100)
            y_smooth = cs(x_smooth)
            ax.plot(x_smooth, y_smooth, "-", color=ax.lines[-1].get_color())

        ax.set_xlabel("Reaction coordinate")
        ax.set_ylabel("Energy (eV)")

        if label is not None:
            ax.legend()

        # Calculate and show energy values in title
        if len(y_valid) >= 2:
            E_barrier = np.max(y_valid)  # Already referenced to initial
            delta_E = y_valid[-1] - y_valid[0]  # Final - initial
            ax.set_title(f"$E_a$ = {E_barrier:.2f} eV, $\\Delta E$ = {delta_E:.2f} eV")

        if show:
            plt.show()

        return ax
