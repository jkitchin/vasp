"""Analysis mixin for post-processing VASP calculations.

Provides methods for:
- Bader charge analysis
- Elastic moduli extraction
- Band structure calculation
- Charge density analysis
"""

from __future__ import annotations

import os
import re
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
        self.atoms.get_cell()
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

    def read_doscar(self) -> dict:
        """Read density of states from DOSCAR.

        Returns:
            Dict with:
            - 'energy': Energy grid in eV
            - 'total_dos': Total DOS
            - 'integrated_dos': Integrated DOS
            - 'fermi': Fermi energy
            - 'projected': Projected DOS per atom (if LORBIT set)

        Raises:
            FileNotFoundError: If DOSCAR not found.
        """
        doscar = os.path.join(self.directory, 'DOSCAR')

        if not os.path.exists(doscar):
            raise FileNotFoundError(f"DOSCAR not found in {self.directory}")

        with open(doscar) as f:
            lines = f.readlines()

        # First line: NIONS, NIONS, JOBPAR, NCDIJ
        # Line 6: EMAX, EMIN, NEDOS, EFERMI, 1.0
        header_parts = lines[5].split()
        emax = float(header_parts[0])
        emin = float(header_parts[1])
        nedos = int(header_parts[2])
        efermi = float(header_parts[3])

        result = {
            'fermi': efermi,
            'emax': emax,
            'emin': emin,
            'nedos': nedos,
        }

        # Parse total DOS (starts at line 6)
        energies = []
        total_dos = []
        integrated_dos = []

        for i in range(nedos):
            line = lines[6 + i]
            parts = line.split()
            energies.append(float(parts[0]))
            total_dos.append(float(parts[1]))
            if len(parts) > 2:
                integrated_dos.append(float(parts[2]))

        result['energy'] = np.array(energies)
        result['total_dos'] = np.array(total_dos)
        if integrated_dos:
            result['integrated_dos'] = np.array(integrated_dos)

        # Check for projected DOS (LORBIT != 0)
        pdos_start = 6 + nedos
        if len(lines) > pdos_start + 1:
            # Parse projected DOS for each atom
            natoms = len(self.atoms)
            projected = []

            for atom_idx in range(natoms):
                atom_start = pdos_start + 1 + atom_idx * (nedos + 1)
                if atom_start >= len(lines):
                    break

                atom_dos = []
                for i in range(nedos):
                    if atom_start + 1 + i >= len(lines):
                        break
                    parts = lines[atom_start + 1 + i].split()
                    # Format depends on LORBIT value
                    # LORBIT=10: s, p, d
                    # LORBIT=11: s, px, py, pz, dxy, dyz, etc.
                    atom_dos.append([float(x) for x in parts[1:]])

                if atom_dos:
                    projected.append(np.array(atom_dos))

            if projected:
                result['projected'] = projected

        return result

    def read_procar(self, kpoint: int = 0, band: int = 0) -> dict:
        """Read projected band character from PROCAR.

        Requires LORBIT=11 or LORBIT=12.

        Args:
            kpoint: K-point index (0-based).
            band: Band index (0-based).

        Returns:
            Dict with:
            - 'kpoints': K-point coordinates
            - 'weights': K-point weights
            - 'energies': Band energies [nkpts, nbands]
            - 'occupations': Occupations [nkpts, nbands]
            - 'projections': Orbital projections [nkpts, nbands, natoms, norbitals]

        Raises:
            FileNotFoundError: If PROCAR not found.
        """
        procar = os.path.join(self.directory, 'PROCAR')

        if not os.path.exists(procar):
            raise FileNotFoundError(f"PROCAR not found in {self.directory}")

        with open(procar) as f:
            content = f.read()

        # Parse header
        header_match = re.search(
            r'# of k-points:\s*(\d+)\s*# of bands:\s*(\d+)\s*# of ions:\s*(\d+)',
            content
        )

        if not header_match:
            raise ValueError("Could not parse PROCAR header")

        nkpts = int(header_match.group(1))
        nbands = int(header_match.group(2))
        nions = int(header_match.group(3))

        result = {
            'nkpts': nkpts,
            'nbands': nbands,
            'nions': nions,
            'kpoints': [],
            'weights': [],
            'energies': np.zeros((nkpts, nbands)),
            'occupations': np.zeros((nkpts, nbands)),
        }

        # Parse k-points
        kpt_pattern = r'k-point\s+(\d+)\s*:\s*([\d.\s-]+)\s+weight\s*=\s*([\d.]+)'
        for match in re.finditer(kpt_pattern, content):
            kpt_coords = [float(x) for x in match.group(2).split()]
            weight = float(match.group(3))
            result['kpoints'].append(kpt_coords)
            result['weights'].append(weight)

        result['kpoints'] = np.array(result['kpoints'])
        result['weights'] = np.array(result['weights'])

        # Parse band energies
        band_pattern = r'band\s+(\d+)\s*#\s*energy\s+([\d.E+-]+)\s*#\s*occ\.\s*([\d.E+-]+)'
        band_matches = list(re.finditer(band_pattern, content))

        for i, match in enumerate(band_matches):
            kpt_idx = i // nbands
            band_idx = i % nbands
            if kpt_idx < nkpts:
                result['energies'][kpt_idx, band_idx] = float(match.group(2))
                result['occupations'][kpt_idx, band_idx] = float(match.group(3))

        return result

    def get_work_function(self, axis: int = 2) -> float:
        """Calculate work function from LOCPOT.

        For a slab calculation with vacuum, calculates the difference
        between vacuum level and Fermi level.

        Args:
            axis: Surface normal direction (0=x, 1=y, 2=z).

        Returns:
            Work function in eV.

        Raises:
            FileNotFoundError: If LOCPOT not found.
        """
        grid, potential = self.get_local_potential()

        # Average potential along surface normal
        axes_to_avg = [i for i in range(3) if i != axis]
        avg_potential = np.mean(potential, axis=tuple(axes_to_avg))

        # Find vacuum level (maximum in vacuum region)
        # Assumes slab is centered
        n_points = len(avg_potential)
        edge_region = n_points // 4

        vacuum_left = np.mean(avg_potential[:edge_region])
        vacuum_right = np.mean(avg_potential[-edge_region:])
        vacuum_level = (vacuum_left + vacuum_right) / 2

        # Get Fermi level from results
        fermi = self.results.get('fermi_level')
        if fermi is None:
            fermi = self._read_fermi_from_vasprun()

        return vacuum_level - fermi

    def get_band_gap_from_doscar(self, tol: float = 0.01) -> float:
        """Calculate band gap from DOSCAR.

        Args:
            tol: DOS threshold for finding band edges.

        Returns:
            Band gap in eV (0 for metals).
        """
        dos_data = self.read_doscar()

        energy = dos_data['energy']
        dos = dos_data['total_dos']
        fermi = dos_data['fermi']

        # Shift to Fermi level
        energy = energy - fermi

        # Find occupied and unoccupied regions
        below_fermi = energy < 0
        above_fermi = energy > 0

        # Find VBM (highest energy below Fermi with non-zero DOS)
        occ_mask = below_fermi & (dos > tol)
        if not np.any(occ_mask):
            return 0.0  # Metal

        vbm_idx = np.where(occ_mask)[0][-1]
        vbm = energy[vbm_idx]

        # Find CBM (lowest energy above Fermi with non-zero DOS)
        unocc_mask = above_fermi & (dos > tol)
        if not np.any(unocc_mask):
            return 0.0  # Metal

        cbm_idx = np.where(unocc_mask)[0][0]
        cbm = energy[cbm_idx]

        gap = cbm - vbm
        return max(0.0, gap)

    def get_elf(self) -> tuple[np.ndarray, np.ndarray]:
        """Read electron localization function from ELFCAR.

        Returns:
            Tuple of (grid_points, elf) arrays.

        Raises:
            FileNotFoundError: If ELFCAR not found.
        """
        elfcar = os.path.join(self.directory, 'ELFCAR')

        if not os.path.exists(elfcar):
            raise FileNotFoundError(f"ELFCAR not found in {self.directory}")

        return self._read_volumetric_file(elfcar)


