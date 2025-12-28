#!/usr/bin/env python
"""
18 - 3D Visualization

This script demonstrates visualization of volumetric data from VASP:
- Charge density (CHGCAR)
- Electrostatic potential (LOCPOT)
- Electron localization function (ELFCAR)

Usage:
    python run.py
"""

import os
import numpy as np
from ase.build import bulk, molecule
from vasp import Vasp

# Try to import matplotlib for plotting
try:
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from mpl_toolkits.mplot3d import Axes3D
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Note: matplotlib not found. Plots will be skipped.")

# Try to import plotly for interactive 3D
try:
    import plotly.graph_objects as go
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False
    print("Note: plotly not found. Interactive 3D plots will be skipped.")


def read_vasp_volumetric(filename):
    """Read VASP volumetric data (CHGCAR, LOCPOT, ELFCAR format)."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Skip header (system name, scale, lattice, atoms)
    # Find the grid dimensions line
    idx = 0
    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) == 3:
            try:
                ngx, ngy, ngz = int(parts[0]), int(parts[1]), int(parts[2])
                idx = i + 1
                break
            except ValueError:
                continue

    # Read data
    data = []
    for line in lines[idx:]:
        # Stop at augmentation charges or spin section
        if 'augmentation' in line.lower():
            break
        values = line.split()
        for v in values:
            try:
                data.append(float(v))
            except ValueError:
                break
        if len(data) >= ngx * ngy * ngz:
            break

    # Reshape to 3D grid
    data = np.array(data[:ngx * ngy * ngz])
    data = data.reshape((ngx, ngy, ngz), order='F')

    return data, (ngx, ngy, ngz)


def plot_2d_slice(data, axis=2, slice_idx=None, title='', filename=None):
    """Plot a 2D slice through volumetric data."""
    if slice_idx is None:
        slice_idx = data.shape[axis] // 2

    if axis == 0:
        slice_data = data[slice_idx, :, :]
        xlabel, ylabel = 'y', 'z'
    elif axis == 1:
        slice_data = data[:, slice_idx, :]
        xlabel, ylabel = 'x', 'z'
    else:
        slice_data = data[:, :, slice_idx]
        xlabel, ylabel = 'x', 'y'

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.contourf(slice_data.T, levels=50, cmap='viridis')
    plt.colorbar(im, ax=ax, label='Value')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_aspect('equal')

    if filename:
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"Saved: {filename}")

    plt.close()


def plot_isosurface_matplotlib(data, isovalue, title='', filename=None):
    """Plot isosurface using matplotlib (simple version)."""
    from skimage import measure

    # Use marching cubes to get isosurface
    try:
        verts, faces, _, _ = measure.marching_cubes(data, isovalue)
    except Exception:
        print(f"Could not generate isosurface at value {isovalue}")
        return

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2],
                    cmap='viridis', alpha=0.7)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)

    if filename:
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"Saved: {filename}")

    plt.close()


def plot_isosurface_plotly(data, isovalue, title='', filename=None):
    """Plot isosurface using plotly (interactive)."""
    X, Y, Z = np.mgrid[0:data.shape[0], 0:data.shape[1], 0:data.shape[2]]

    fig = go.Figure(data=go.Isosurface(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=data.flatten(),
        isomin=isovalue * 0.9,
        isomax=isovalue * 1.1,
        surface_count=2,
        colorscale='Viridis',
        caps=dict(x_show=False, y_show=False, z_show=False),
        opacity=0.6,
    ))

    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            aspectmode='data'
        )
    )

    if filename:
        fig.write_html(filename)
        print(f"Saved: {filename}")

    return fig


# =============================================================================
# Main Calculation
# =============================================================================

print("=" * 60)
print("3D Volumetric Data Visualization")
print("=" * 60)
print()

# =============================================================================
# Step 1: Run VASP calculation with volumetric output
# =============================================================================

print("Step 1: Running VASP calculation")
print()

# Use silicon as example
atoms = bulk('Si', 'diamond', a=5.43)

calc = Vasp(
    label='viz/si',
    atoms=atoms,
    xc='PBE',
    encut=300,      # Lower for speed
    kpts=(4, 4, 4),
    ismear=0,
    sigma=0.1,
    lcharg=True,    # Write CHGCAR
    lvtot=True,     # Write LOCPOT
    lelf=True,      # Write ELFCAR
    prec='Normal',
)

energy = calc.potential_energy
print(f"  Total energy: {energy:.6f} eV")
print()

# =============================================================================
# Step 2: Read volumetric data
# =============================================================================

print("Step 2: Reading volumetric data")
print()

calc_dir = calc.directory

# Check which files exist
files_to_read = ['CHGCAR', 'LOCPOT', 'ELFCAR']
available_files = {}

for fname in files_to_read:
    fpath = os.path.join(calc_dir, fname)
    if os.path.exists(fpath):
        try:
            data, grid = read_vasp_volumetric(fpath)
            available_files[fname] = data
            print(f"  {fname}: grid = {grid}, range = [{data.min():.3f}, {data.max():.3f}]")
        except Exception as e:
            print(f"  {fname}: could not read ({e})")
    else:
        print(f"  {fname}: not found (file will be generated by VASP)")

print()

# =============================================================================
# Step 3: Visualize charge density
# =============================================================================

if HAS_MATPLOTLIB:
    print("Step 3: Visualizing volumetric data")
    print()

    if 'CHGCAR' in available_files:
        chg = available_files['CHGCAR']

        # Normalize by volume
        vol = atoms.cell.volume
        chg_normalized = chg / vol

        # 2D slices
        plot_2d_slice(
            chg_normalized,
            axis=2,
            title='Charge Density - XY plane (z=0.5)',
            filename='charge_density_xy.png'
        )

        plot_2d_slice(
            chg_normalized,
            axis=0,
            title='Charge Density - YZ plane (x=0.5)',
            filename='charge_density_yz.png'
        )

        # Find appropriate isosurface value
        mean_chg = chg_normalized.mean()
        isovalue = mean_chg * 2  # Above average density

        print(f"  Charge density: mean = {mean_chg:.4f} e/Å³")
        print(f"  Using isosurface value: {isovalue:.4f} e/Å³")

        # Try isosurface with scikit-image
        try:
            from skimage import measure
            plot_isosurface_matplotlib(
                chg_normalized,
                isovalue,
                title=f'Si Charge Density Isosurface (ρ = {isovalue:.3f} e/Å³)',
                filename='charge_density_3d.png'
            )
        except ImportError:
            print("  Note: scikit-image not found, skipping 3D isosurface")

    if 'LOCPOT' in available_files:
        pot = available_files['LOCPOT']

        plot_2d_slice(
            pot,
            axis=2,
            title='Local Potential - XY plane (z=0.5)',
            filename='potential_xy.png'
        )

        print(f"  Local potential: range = [{pot.min():.2f}, {pot.max():.2f}] eV")

    if 'ELFCAR' in available_files:
        elf = available_files['ELFCAR']

        plot_2d_slice(
            elf,
            axis=2,
            title='ELF - XY plane (z=0.5)',
            filename='elf_xy.png'
        )

        # ELF isosurface at high localization
        try:
            from skimage import measure
            plot_isosurface_matplotlib(
                elf,
                isovalue=0.8,
                title='Si ELF Isosurface (ELF = 0.8)',
                filename='elf_3d.png'
            )
        except ImportError:
            pass

        print(f"  ELF: range = [{elf.min():.3f}, {elf.max():.3f}]")

    print()

# =============================================================================
# Step 4: Interactive 3D with Plotly (if available)
# =============================================================================

if HAS_PLOTLY and 'CHGCAR' in available_files:
    print("Step 4: Creating interactive 3D visualization")
    print()

    chg = available_files['CHGCAR']
    vol = atoms.cell.volume
    chg_normalized = chg / vol
    isovalue = chg_normalized.mean() * 2

    fig = plot_isosurface_plotly(
        chg_normalized,
        isovalue,
        title=f'Si Charge Density (interactive)',
        filename='charge_density_interactive.html'
    )

    print()

# =============================================================================
# Step 5: Generate mock data for demonstration (if no VASP output)
# =============================================================================

if not available_files and HAS_MATPLOTLIB:
    print("Step 5: Generating demonstration with mock data")
    print()

    # Create mock charge density (Gaussian at atom positions)
    ngx, ngy, ngz = 50, 50, 50
    mock_chg = np.zeros((ngx, ngy, ngz))

    # Add Gaussian blobs at approximate Si positions in diamond
    positions = [
        (0, 0, 0), (0.25, 0.25, 0.25),
        (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5),
        (0.75, 0.75, 0.25), (0.75, 0.25, 0.75), (0.25, 0.75, 0.75)
    ]

    for pos in positions:
        x0, y0, z0 = int(pos[0] * ngx), int(pos[1] * ngy), int(pos[2] * ngz)
        for i in range(ngx):
            for j in range(ngy):
                for k in range(ngz):
                    # Distance with periodic boundary
                    dx = min(abs(i - x0), ngx - abs(i - x0))
                    dy = min(abs(j - y0), ngy - abs(j - y0))
                    dz = min(abs(k - z0), ngz - abs(k - z0))
                    r2 = dx**2 + dy**2 + dz**2
                    mock_chg[i, j, k] += np.exp(-r2 / 20)

    plot_2d_slice(
        mock_chg,
        axis=2,
        title='Mock Charge Density - XY plane',
        filename='mock_charge_xy.png'
    )

    print("  Generated mock visualization")
    print()

# =============================================================================
# Summary
# =============================================================================

print("=" * 60)
print("3D Visualization complete!")
print("=" * 60)
print()
print("Generated files:")

for ext in ['png', 'html']:
    import glob
    files = glob.glob(f'*.{ext}')
    for f in files:
        print(f"  - {f}")

print()
print("Key points:")
print("  - CHGCAR: electron density, analyze bonding")
print("  - LOCPOT: electrostatic potential, work function")
print("  - ELFCAR: electron localization, lone pairs/bonds")
print("  - Use VESTA or ParaView for publication-quality rendering")
print()
print("Recommended tools for better visualization:")
print("  - VESTA: https://jp-minerals.org/vesta/")
print("  - ParaView: https://www.paraview.org/")
print("  - py3Dmol: pip install py3Dmol")
print()
print("Congratulations! You've completed the tutorial series.")
