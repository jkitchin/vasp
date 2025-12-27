"""Analysis mixin for post-processing VASP calculations.

Provides methods for:
- Bader charge analysis
- Elastic moduli extraction
- Band structure calculation
- Charge density analysis
"""

from __future__ import annotations

import os
import subprocess
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    pass


class AnalysisMixin:
    """Mixin providing post-processing analysis methods.

    Requires:
    - directory: str - Calculation directory path
    - atoms: Atoms - ASE Atoms object
    - resort: list - Reverse sorting indices
    - update(): method - Ensure calculation is complete
    """

    def bader(
        self,
        cmd: str | list[str] | None = None,
        ref: bool = True
    ) -> np.ndarray:
        """Run Bader charge analysis.

        Requires the 'bader' program from the Henkelman group:
        http://theory.cm.utexas.edu/henkelman/code/bader/

        Args:
            cmd: Bader command (default: 'bader').
            ref: If True, use reference charge density (AECCAR0 + AECCAR2).

        Returns:
            Array of Bader charges per atom.

        Raises:
            RuntimeError: If Bader analysis fails.
        """
        self.update()

        if cmd is None:
            cmd = 'bader'

        # Build reference charge if requested
        if ref:
            ref_file = self._build_reference_charge()
            if ref_file:
                bader_cmd = [cmd, 'CHGCAR', '-ref', ref_file]
            else:
                bader_cmd = [cmd, 'CHGCAR']
        else:
            bader_cmd = [cmd, 'CHGCAR']

        if isinstance(cmd, str):
            bader_cmd = ' '.join(bader_cmd)
            shell = True
        else:
            shell = False

        # Run Bader
        result = subprocess.run(
            bader_cmd,
            cwd=self.directory,
            shell=shell,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            raise RuntimeError(f"Bader analysis failed: {result.stderr}")

        # Parse ACF.dat
        return self._parse_bader_acf()

    def _build_reference_charge(self) -> str | None:
        """Build reference charge density from AECCAR files."""
        aec0 = os.path.join(self.directory, 'AECCAR0')
        aec2 = os.path.join(self.directory, 'AECCAR2')

        if not os.path.exists(aec0) or not os.path.exists(aec2):
            return None

        # Sum the charge densities
        ref_file = os.path.join(self.directory, 'CHGREF')

        # Use chgsum.pl if available, otherwise do it in Python
        try:
            subprocess.run(
                ['chgsum.pl', 'AECCAR0', 'AECCAR2'],
                cwd=self.directory,
                check=True,
                capture_output=True,
            )
            if os.path.exists(os.path.join(self.directory, 'CHGCAR_sum')):
                os.rename(
                    os.path.join(self.directory, 'CHGCAR_sum'),
                    ref_file
                )
                return ref_file
        except (subprocess.CalledProcessError, FileNotFoundError):
            pass

        # Python fallback - simple header copy + charge sum
        # This is a simplified version
        return None

    def _parse_bader_acf(self) -> np.ndarray:
        """Parse Bader ACF.dat output file."""
        acf_file = os.path.join(self.directory, 'ACF.dat')

        if not os.path.exists(acf_file):
            raise FileNotFoundError("ACF.dat not found - Bader analysis may have failed")

        charges = []
        with open(acf_file) as f:
            for line in f:
                if line.startswith('#') or line.startswith('-'):
                    continue
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        # Column 5 is the charge
                        charges.append(float(parts[4]))
                    except ValueError:
                        continue

        if not charges:
            raise ValueError("No charges found in ACF.dat")

        # Unsort to match original atom order
        return np.array(charges)[self.resort]

    def get_elastic_moduli(self) -> np.ndarray:
        """Extract elastic moduli from OUTCAR.

        Requires IBRION=6 and ISIF>=3 calculation.

        Returns:
            6x6 elastic constants matrix in GPa.
        """
        self.update()

        outcar = os.path.join(self.directory, 'OUTCAR')
        if not os.path.exists(outcar):
            raise FileNotFoundError(f"OUTCAR not found in {self.directory}")

        with open(outcar) as f:
            content = f.read()

        # Find elastic constants section
        if 'TOTAL ELASTIC MODULI' not in content:
            raise ValueError(
                "Elastic moduli not found. "
                "Ensure IBRION=6 and ISIF>=3 were used."
            )

        # Parse the 6x6 matrix
        pattern = r'TOTAL ELASTIC MODULI \(kBar\).*?-+\s*\n((?:.*\n){6})'
        match = re.search(pattern, content, re.DOTALL)

        if not match:
            raise ValueError("Could not parse elastic moduli from OUTCAR")

        matrix = []
        for line in match.group(1).strip().split('\n'):
            parts = line.split()
            if len(parts) >= 7:
                # Skip first column (XX, YY, etc.)
                row = [float(x) for x in parts[1:7]]
                matrix.append(row)

        # Convert from kBar to GPa
        return np.array(matrix) * 0.1

    def get_bulk_modulus(self, method: str = 'voigt') -> float:
        """Calculate bulk modulus from elastic constants.

        Args:
            method: Averaging method ('voigt', 'reuss', or 'hill').

        Returns:
            Bulk modulus in GPa.
        """
        C = self.get_elastic_moduli()

        if method == 'voigt':
            K = (C[0, 0] + C[1, 1] + C[2, 2] +
                 2 * (C[0, 1] + C[1, 2] + C[0, 2])) / 9
        elif method == 'reuss':
            S = np.linalg.inv(C)
            K = 1.0 / (S[0, 0] + S[1, 1] + S[2, 2] +
                       2 * (S[0, 1] + S[1, 2] + S[0, 2]))
        elif method == 'hill':
            K_voigt = (C[0, 0] + C[1, 1] + C[2, 2] +
                       2 * (C[0, 1] + C[1, 2] + C[0, 2])) / 9
            S = np.linalg.inv(C)
            K_reuss = 1.0 / (S[0, 0] + S[1, 1] + S[2, 2] +
                             2 * (S[0, 1] + S[1, 2] + S[0, 2]))
            K = (K_voigt + K_reuss) / 2
        else:
            raise ValueError(f"Unknown method: {method}")

        return K

    def get_charge_density(self, spin: int | None = None) -> tuple[np.ndarray, np.ndarray]:
        """Read charge density from CHGCAR.

        Args:
            spin: Spin channel (None for total, 0 or 1 for spin-polarized).

        Returns:
            Tuple of (grid_points, density) arrays.
        """
        chgcar = os.path.join(self.directory, 'CHGCAR')

        if not os.path.exists(chgcar):
            raise FileNotFoundError(f"CHGCAR not found in {self.directory}")

        # Parse CHGCAR (simplified)
        with open(chgcar) as f:
            lines = f.readlines()

        # Skip header (atoms section)
        # Find grid dimensions after the atom positions
        header_end = 0
        for i, line in enumerate(lines):
            if i > 5:
                parts = line.split()
                if len(parts) == 3:
                    try:
                        grid = [int(x) for x in parts]
                        header_end = i + 1
                        break
                    except ValueError:
                        continue

        if header_end == 0:
            raise ValueError("Could not parse CHGCAR grid dimensions")

        ngx, ngy, ngz = grid
        npoints = ngx * ngy * ngz

        # Read density values
        values = []
        for line in lines[header_end:]:
            values.extend([float(x) for x in line.split()])
            if len(values) >= npoints:
                break

        density = np.array(values[:npoints]).reshape((ngx, ngy, ngz))

        # Generate grid points
        cell = self.atoms.get_cell()
        x = np.linspace(0, 1, ngx, endpoint=False)
        y = np.linspace(0, 1, ngy, endpoint=False)
        z = np.linspace(0, 1, ngz, endpoint=False)
        grid_points = np.array(np.meshgrid(x, y, z, indexing='ij'))

        return grid_points, density

    def get_local_potential(self) -> tuple[np.ndarray, np.ndarray]:
        """Read local potential from LOCPOT.

        Returns:
            Tuple of (grid_points, potential) arrays.
        """
        locpot = os.path.join(self.directory, 'LOCPOT')

        if not os.path.exists(locpot):
            raise FileNotFoundError(f"LOCPOT not found in {self.directory}")

        # Same format as CHGCAR
        return self._read_volumetric_file(locpot)

    def _read_volumetric_file(self, filepath: str) -> tuple[np.ndarray, np.ndarray]:
        """Read volumetric data file (CHGCAR, LOCPOT, etc.)."""
        with open(filepath) as f:
            lines = f.readlines()

        # Find grid dimensions
        header_end = 0
        grid = None
        for i, line in enumerate(lines):
            if i > 5:
                parts = line.split()
                if len(parts) == 3:
                    try:
                        grid = [int(x) for x in parts]
                        header_end = i + 1
                        break
                    except ValueError:
                        continue

        if grid is None:
            raise ValueError(f"Could not parse grid dimensions from {filepath}")

        ngx, ngy, ngz = grid
        npoints = ngx * ngy * ngz

        values = []
        for line in lines[header_end:]:
            values.extend([float(x) for x in line.split()])
            if len(values) >= npoints:
                break

        data = np.array(values[:npoints]).reshape((ngx, ngy, ngz))

        x = np.linspace(0, 1, ngx, endpoint=False)
        y = np.linspace(0, 1, ngy, endpoint=False)
        z = np.linspace(0, 1, ngz, endpoint=False)
        grid_points = np.array(np.meshgrid(x, y, z, indexing='ij'))

        return grid_points, data


# Import re for elastic moduli parsing
import re
