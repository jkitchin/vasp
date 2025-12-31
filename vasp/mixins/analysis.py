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
        self, cmd: str | list[str] | None = None, ref: bool = True, overwrite: bool = False
    ) -> np.ndarray:
        """Run Bader charge analysis.

        Requires the 'bader' program from the Henkelman group:
        http://theory.cm.utexas.edu/henkelman/code/bader/

        Args:
            cmd: Bader command (default: 'bader').
            ref: If True, use reference charge density (AECCAR0 + AECCAR2).
            overwrite: If False and ACF.dat exists, use existing results.
                If True, always run bader analysis.

        Returns:
            Array of Bader charges per atom.

        Raises:
            RuntimeError: If Bader analysis fails.
        """
        acf_file = os.path.join(self.directory, "ACF.dat")

        # Skip running bader if ACF.dat exists and overwrite=False
        if not overwrite and os.path.exists(acf_file):
            pass  # Just parse existing results
        else:
            # Need to run bader - ensure calculation is complete first
            self.update()
            if cmd is None:
                cmd = "bader"

            # Build reference charge if requested
            if ref:
                ref_file = self._build_reference_charge()
                if ref_file:
                    bader_cmd = [cmd, "CHGCAR", "-ref", ref_file]
                else:
                    bader_cmd = [cmd, "CHGCAR"]
            else:
                bader_cmd = [cmd, "CHGCAR"]

            if isinstance(cmd, str):
                bader_cmd = " ".join(bader_cmd)
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

        # Parse ACF.dat - raw Bader electron counts
        bader_electrons = self._parse_bader_acf()

        # Convert to ionic charges: charge = ZVAL - bader_electrons
        zval_per_atom = self._get_zval_per_atom()
        charges = zval_per_atom - bader_electrons

        # Store charges in results for load_atoms() to apply
        self.results["bader_charges"] = charges

        # Set charges on atoms if already loaded
        if self.atoms is not None:
            self.atoms.set_initial_charges(charges)

        return charges

    def _build_reference_charge(self) -> str | None:
        """Build reference charge density from AECCAR files."""
        aec0 = os.path.join(self.directory, "AECCAR0")
        aec2 = os.path.join(self.directory, "AECCAR2")

        if not os.path.exists(aec0) or not os.path.exists(aec2):
            return None

        ref_file = os.path.join(self.directory, "CHGREF")

        # If CHGREF already exists, use it
        if os.path.exists(ref_file):
            return ref_file

        # Use chgsum.pl if available, otherwise do it in Python
        try:
            subprocess.run(
                ["chgsum.pl", "AECCAR0", "AECCAR2"],
                cwd=self.directory,
                check=True,
                capture_output=True,
            )
            if os.path.exists(os.path.join(self.directory, "CHGCAR_sum")):
                os.rename(os.path.join(self.directory, "CHGCAR_sum"), ref_file)
                return ref_file
        except (subprocess.CalledProcessError, FileNotFoundError):
            pass

        # Python fallback - sum AECCAR0 + AECCAR2
        try:
            self._sum_charge_files(aec0, aec2, ref_file)
            return ref_file
        except Exception:
            return None

    def _sum_charge_files(self, file1: str, file2: str, output: str) -> None:
        """Sum two VASP charge density files (AECCAR0 + AECCAR2)."""
        with open(file1) as f:
            lines1 = f.readlines()
        with open(file2) as f:
            lines2 = f.readlines()

        # Find where the charge data starts (after blank line following atom positions)
        header_end = 0
        for i, line in enumerate(lines1):
            if i > 5:
                stripped = line.strip()
                if stripped == "" and i + 1 < len(lines1):
                    next_parts = lines1[i + 1].split()
                    if len(next_parts) == 3:
                        try:
                            grid = [int(x) for x in next_parts]
                            if all(g > 1 for g in grid):
                                header_end = i + 2
                                break
                        except ValueError:
                            continue

        # Write header unchanged
        with open(output, "w") as out:
            for line in lines1[:header_end]:
                out.write(line)

            # Sum the charge data
            for l1, l2 in zip(lines1[header_end:], lines2[header_end:]):
                vals1 = l1.split()
                vals2 = l2.split()
                if len(vals1) == len(vals2) and vals1:
                    try:
                        summed = [float(v1) + float(v2) for v1, v2 in zip(vals1, vals2)]
                        out.write(" ".join(f"{v:.11E}" for v in summed) + "\n")
                    except ValueError:
                        # Not numeric data, write as-is
                        out.write(l1)
                else:
                    out.write(l1)

    def _parse_bader_acf(self) -> np.ndarray:
        """Parse Bader ACF.dat output file."""
        acf_file = os.path.join(self.directory, "ACF.dat")

        if not os.path.exists(acf_file):
            raise FileNotFoundError("ACF.dat not found - Bader analysis may have failed")

        charges = []
        with open(acf_file) as f:
            for line in f:
                if line.startswith("#") or line.startswith("-"):
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

        charges_array = np.array(charges)

        # Unsort to match original atom order (if resort is available)
        if self.resort:
            return charges_array[self.resort]
        return charges_array

    def _get_zval_per_atom(self) -> np.ndarray:
        """Get valence electrons (ZVAL) for each atom from POTCAR.

        Returns:
            Array of ZVAL values, one per atom in the structure.
        """
        from ase.io import read

        potcar = os.path.join(self.directory, "POTCAR")
        if not os.path.exists(potcar):
            raise FileNotFoundError(f"POTCAR not found in {self.directory}")

        # Read ZVAL for each species from POTCAR
        zvals = []
        with open(potcar) as f:
            for line in f:
                if "ZVAL" in line:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == "ZVAL":
                            zvals.append(float(parts[i + 2]))
                            break

        # Read atoms to get the species order and counts
        contcar = os.path.join(self.directory, "CONTCAR")
        poscar = os.path.join(self.directory, "POSCAR")

        if os.path.exists(contcar) and os.path.getsize(contcar) > 0:
            atoms = read(contcar, format="vasp")
        elif os.path.exists(poscar):
            atoms = read(poscar, format="vasp")
        else:
            raise FileNotFoundError(f"No POSCAR/CONTCAR in {self.directory}")

        # Get species order from structure
        symbols = atoms.get_chemical_symbols()
        species_order = []
        for s in symbols:
            if s not in species_order:
                species_order.append(s)

        # Build ZVAL per atom array
        zval_per_atom = []
        zval_dict = dict(zip(species_order, zvals))
        for symbol in symbols:
            zval_per_atom.append(zval_dict[symbol])

        return np.array(zval_per_atom)

    def get_elastic_moduli(self) -> np.ndarray:
        """Extract elastic moduli from OUTCAR.

        Requires IBRION=6 and ISIF>=3 calculation.

        Returns:
            6x6 elastic constants matrix in GPa.
        """
        self.update()

        outcar = os.path.join(self.directory, "OUTCAR")
        if not os.path.exists(outcar):
            raise FileNotFoundError(f"OUTCAR not found in {self.directory}")

        with open(outcar) as f:
            lines = f.readlines()

        # Find TOTAL ELASTIC MODULI section (more robust line-by-line approach)
        matrix = []
        for i, line in enumerate(lines):
            if "TOTAL ELASTIC MODULI (kBar)" in line:
                # Skip header line, "Direction XX YY..." line, and dashed line
                # The matrix data starts 3 lines after "TOTAL ELASTIC MODULI"
                data_lines = lines[i + 3 : i + 9]  # 6 rows of data
                for data_line in data_lines:
                    parts = data_line.split()
                    # Each line: XX/YY/ZZ/XY/YZ/ZX followed by 6 numbers
                    if len(parts) >= 7 and parts[0] in ("XX", "YY", "ZZ", "XY", "YZ", "ZX"):
                        row = [float(x) for x in parts[1:7]]
                        matrix.append(row)
                break

        if len(matrix) != 6:
            raise ValueError(
                "Elastic moduli not found or incomplete. " "Ensure IBRION=6 and ISIF>=3 were used."
            )

        # Convert from kBar to GPa
        return np.array(matrix) * 0.1

    def get_bulk_modulus(self, method: str = "voigt") -> float:
        """Calculate bulk modulus from elastic constants.

        Args:
            method: Averaging method ('voigt', 'reuss', or 'hill').

        Returns:
            Bulk modulus in GPa.
        """
        C = self.get_elastic_moduli()

        if method == "voigt":
            K = (C[0, 0] + C[1, 1] + C[2, 2] + 2 * (C[0, 1] + C[1, 2] + C[0, 2])) / 9
        elif method == "reuss":
            S = np.linalg.inv(C)
            K = 1.0 / (S[0, 0] + S[1, 1] + S[2, 2] + 2 * (S[0, 1] + S[1, 2] + S[0, 2]))
        elif method == "hill":
            K_voigt = (C[0, 0] + C[1, 1] + C[2, 2] + 2 * (C[0, 1] + C[1, 2] + C[0, 2])) / 9
            S = np.linalg.inv(C)
            K_reuss = 1.0 / (S[0, 0] + S[1, 1] + S[2, 2] + 2 * (S[0, 1] + S[1, 2] + S[0, 2]))
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
        chgcar = os.path.join(self.directory, "CHGCAR")

        if not os.path.exists(chgcar):
            raise FileNotFoundError(f"CHGCAR not found in {self.directory}")

        if os.path.getsize(chgcar) == 0:
            raise ValueError(
                "CHGCAR is empty. Set lcharg=True in your calculation to write charge density."
            )

        # Use shared volumetric file reader
        return self._read_volumetric_file(chgcar)

    def get_local_potential(self) -> tuple[np.ndarray, np.ndarray]:
        """Read local potential from LOCPOT.

        Returns:
            Tuple of (grid_points, potential) arrays.
        """
        locpot = os.path.join(self.directory, "LOCPOT")

        if not os.path.exists(locpot):
            raise FileNotFoundError(f"LOCPOT not found in {self.directory}")

        # Same format as CHGCAR
        return self._read_volumetric_file(locpot)

    def _read_volumetric_file(self, filepath: str) -> tuple[np.ndarray, np.ndarray]:
        """Read volumetric data file (CHGCAR, LOCPOT, etc.)."""
        with open(filepath) as f:
            lines = f.readlines()

        # Find grid dimensions - look for blank line followed by 3 integers
        # The format is: POSCAR header, blank line, NGX NGY NGZ, then data
        header_end = 0
        grid = None
        for i, line in enumerate(lines):
            if i > 5:
                stripped = line.strip()
                # Look for blank line followed by grid dimensions
                if stripped == "" and i + 1 < len(lines):
                    next_parts = lines[i + 1].split()
                    if len(next_parts) == 3:
                        try:
                            grid = [int(x) for x in next_parts]
                            # Grid dimensions should be reasonably large (> 1)
                            if all(g > 1 for g in grid):
                                header_end = i + 2
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
        grid_points = np.array(np.meshgrid(x, y, z, indexing="ij"))

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
        doscar = os.path.join(self.directory, "DOSCAR")

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
            "fermi": efermi,
            "emax": emax,
            "emin": emin,
            "nedos": nedos,
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

        result["energy"] = np.array(energies)
        result["total_dos"] = np.array(total_dos)
        if integrated_dos:
            result["integrated_dos"] = np.array(integrated_dos)

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
                result["projected"] = projected

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
        procar = os.path.join(self.directory, "PROCAR")

        if not os.path.exists(procar):
            raise FileNotFoundError(f"PROCAR not found in {self.directory}")

        with open(procar) as f:
            content = f.read()

        # Parse header
        header_match = re.search(
            r"# of k-points:\s*(\d+)\s*# of bands:\s*(\d+)\s*# of ions:\s*(\d+)", content
        )

        if not header_match:
            raise ValueError("Could not parse PROCAR header")

        nkpts = int(header_match.group(1))
        nbands = int(header_match.group(2))
        nions = int(header_match.group(3))

        result = {
            "nkpts": nkpts,
            "nbands": nbands,
            "nions": nions,
            "kpoints": [],
            "weights": [],
            "energies": np.zeros((nkpts, nbands)),
            "occupations": np.zeros((nkpts, nbands)),
        }

        # Parse k-points
        kpt_pattern = r"k-point\s+(\d+)\s*:\s*([\d.\s-]+)\s+weight\s*=\s*([\d.]+)"
        for match in re.finditer(kpt_pattern, content):
            kpt_coords = [float(x) for x in match.group(2).split()]
            weight = float(match.group(3))
            result["kpoints"].append(kpt_coords)
            result["weights"].append(weight)

        result["kpoints"] = np.array(result["kpoints"])
        result["weights"] = np.array(result["weights"])

        # Parse band energies
        band_pattern = r"band\s+(\d+)\s*#\s*energy\s+([\d.E+-]+)\s*#\s*occ\.\s*([\d.E+-]+)"
        band_matches = list(re.finditer(band_pattern, content))

        for i, match in enumerate(band_matches):
            kpt_idx = i // nbands
            band_idx = i % nbands
            if kpt_idx < nkpts:
                result["energies"][kpt_idx, band_idx] = float(match.group(2))
                result["occupations"][kpt_idx, band_idx] = float(match.group(3))

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
        fermi = self.results.get("fermi_level")
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

        energy = dos_data["energy"]
        dos = dos_data["total_dos"]
        fermi = dos_data["fermi"]

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
        elfcar = os.path.join(self.directory, "ELFCAR")

        if not os.path.exists(elfcar):
            raise FileNotFoundError(f"ELFCAR not found in {self.directory}")

        return self._read_volumetric_file(elfcar)

    def get_elapsed_time(self) -> float:
        """Get elapsed wall-clock time from OUTCAR.

        Returns:
            Elapsed time in seconds.

        Raises:
            FileNotFoundError: If OUTCAR not found.
            ValueError: If timing information not found.
        """
        import re

        outcar = os.path.join(self.directory, "OUTCAR")

        if not os.path.exists(outcar):
            raise FileNotFoundError(f"OUTCAR not found in {self.directory}")

        with open(outcar) as f:
            content = f.read()

        match = re.search(r"Elapsed time \(sec\):\s*([\d.]+)", content)
        if match:
            return float(match.group(1))

        raise ValueError(f"Elapsed time not found in {outcar}")

    def plot_isosurface(
        self,
        data: np.ndarray | str = "charge",
        isovalue: float | None = None,
        color: str = "blue",
        alpha: float = 0.3,
        show_cell: bool = True,
        show_atoms: bool = True,
        ax=None,
    ):
        """Plot isosurface of 3D volumetric data.

        Args:
            data: 3D numpy array or string ('charge', 'potential', 'elf').
            isovalue: Isosurface value. If None, uses 10% of max value.
            color: Surface color.
            alpha: Surface transparency (0-1).
            show_cell: If True, draw unit cell edges.
            show_atoms: If True, show atom positions.
            ax: Matplotlib 3D axes. If None, creates new figure.

        Returns:
            Matplotlib 3D axes object.

        Raises:
            ImportError: If scikit-image is not installed.
        """
        try:
            from skimage import measure
        except ImportError:
            raise ImportError(
                "scikit-image is required for isosurface plotting. "
                "Install with: pip install scikit-image"
            ) from None

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        # Get the data array
        if isinstance(data, str):
            if data == "charge":
                _, data_array = self.get_charge_density()
            elif data == "potential":
                _, data_array = self.get_local_potential()
            elif data == "elf":
                _, data_array = self.get_elf()
            else:
                raise ValueError(f"Unknown data type: {data}. Use 'charge', 'potential', or 'elf'.")
        else:
            data_array = data

        # Default isovalue
        if isovalue is None:
            isovalue = 0.1 * np.max(data_array)

        # Get unit cell
        from ase.io import read

        contcar = os.path.join(self.directory, "CONTCAR")
        poscar = os.path.join(self.directory, "POSCAR")
        if os.path.exists(contcar) and os.path.getsize(contcar) > 0:
            atoms = read(contcar, format="vasp")
        elif os.path.exists(poscar):
            atoms = read(poscar, format="vasp")
        else:
            atoms = self.atoms

        cell = atoms.get_cell()

        # Extract isosurface using marching cubes
        try:
            verts, faces, normals, values = measure.marching_cubes(
                data_array,
                level=isovalue,
                spacing=(
                    1.0 / data_array.shape[0],
                    1.0 / data_array.shape[1],
                    1.0 / data_array.shape[2],
                ),
            )
        except ValueError:
            raise ValueError(
                f"No isosurface found at isovalue={isovalue}. "
                f"Data range: [{np.min(data_array):.3e}, {np.max(data_array):.3e}]"
            ) from None

        # Transform vertices from fractional to Cartesian coordinates
        verts_cartesian = verts @ cell

        # Create figure if needed
        if ax is None:
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection="3d")

        # Create mesh
        mesh = Poly3DCollection(verts_cartesian[faces], alpha=alpha)
        mesh.set_facecolor(color)
        mesh.set_edgecolor("none")
        ax.add_collection3d(mesh)

        # Draw unit cell
        if show_cell:
            self._draw_cell(ax, cell)

        # Draw atoms
        if show_atoms and atoms is not None:
            positions = atoms.get_positions()
            symbols = atoms.get_chemical_symbols()

            # Simple color mapping
            colors = {
                "H": "white",
                "C": "gray",
                "N": "blue",
                "O": "red",
                "S": "yellow",
                "Fe": "orange",
                "Ni": "green",
                "Cu": "brown",
                "Pt": "silver",
                "Au": "gold",
            }

            for pos, sym in zip(positions, symbols):
                c = colors.get(sym, "purple")
                ax.scatter(*pos, c=c, s=100, edgecolors="black", linewidths=0.5)

        # Set axis limits
        all_verts = np.vstack(
            [
                [0, 0, 0],
                cell[0],
                cell[1],
                cell[2],
                cell[0] + cell[1],
                cell[1] + cell[2],
                cell[0] + cell[2],
                cell[0] + cell[1] + cell[2],
            ]
        )
        margin = 0.5
        ax.set_xlim(all_verts[:, 0].min() - margin, all_verts[:, 0].max() + margin)
        ax.set_ylim(all_verts[:, 1].min() - margin, all_verts[:, 1].max() + margin)
        ax.set_zlim(all_verts[:, 2].min() - margin, all_verts[:, 2].max() + margin)

        ax.set_xlabel("x (Å)")
        ax.set_ylabel("y (Å)")
        ax.set_zlabel("z (Å)")

        return ax

    def _draw_cell(self, ax, cell):
        """Draw unit cell edges on 3D axes."""
        # Cell vertices
        o = np.array([0, 0, 0])
        a, b, c = cell[0], cell[1], cell[2]

        # 12 edges of the parallelepiped
        edges = [
            (o, a),
            (o, b),
            (o, c),
            (a, a + b),
            (a, a + c),
            (b, b + a),
            (b, b + c),
            (c, c + a),
            (c, c + b),
            (a + b, a + b + c),
            (a + c, a + b + c),
            (b + c, a + b + c),
        ]

        for start, end in edges:
            ax.plot3D(*zip(start, end), "k-", linewidth=0.5)
