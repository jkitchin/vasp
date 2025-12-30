"""Command-line interface for vasp-ase.

Provides utilities for:
- Checking VASP job status
- Installing Claude Code skills globally
- Package information
"""

import argparse
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Optional

__all__ = ["main", "status_command", "claude_install", "claude_uninstall", "vaspsum", "vasp_debug"]


def get_claude_home() -> Path:
    """Get the Claude Code configuration directory."""
    # Check for custom location
    claude_home = os.environ.get("CLAUDE_CONFIG_DIR")
    if claude_home:
        return Path(claude_home)

    # Default to ~/.claude
    return Path.home() / ".claude"


def get_package_claude_dir() -> Optional[Path]:
    """Get the .claude directory from the installed package."""
    # Try to find the package's .claude directory
    try:
        import vasp

        package_dir = Path(vasp.__file__).parent.parent
        claude_dir = package_dir / ".claude"
        if claude_dir.exists():
            return claude_dir
    except Exception:
        pass

    # Fallback: try relative to this file
    this_dir = Path(__file__).parent.parent
    claude_dir = this_dir / ".claude"
    if claude_dir.exists():
        return claude_dir

    return None


# Global skill content (embedded for reliable installation)
GLOBAL_SKILL_CONTENT = """---
name: vasp
description: Help with VASP DFT calculations using the vasp-ase Python interface. Use when the user asks about VASP, DFT calculations, or computational materials science.
---

# VASP-ASE Interface Skill

You have access to the vasp-ase package for running VASP calculations.

## Installation Check

First, verify the package is available:
```python
import vasp
print(f"vasp-ase version: {vasp.__version__}")
```

## Quick Start

```python
from ase.build import bulk
from vasp import Vasp

# Create structure
atoms = bulk('Si', 'diamond', a=5.43)

# Run calculation
calc = Vasp(
    atoms=atoms,
    xc='PBE',
    encut=400,
    kpts=(8, 8, 8),
)

energy = calc.potential_energy
```

## Parameter Presets

```python
from vasp.parameters import get_vdw_params, get_ldau_params, get_hybrid_params, HubbardU

# Van der Waals (D3-BJ)
vdw = get_vdw_params('d3bj')  # Also: 'd2', 'd3', 'ts'

# DFT+U
ldau = get_ldau_params(['Fe', 'O'], {'Fe': HubbardU(u=4.0)})

# Hybrid functionals
hse = get_hybrid_params('hse06')  # Also: 'pbe0', 'b3lyp'

# Use with calculator
calc = Vasp(atoms=atoms, **vdw)
```

## Recipes for Common Workflows

```python
from vasp.recipes import static_job, relax_job, double_relax_flow

# Single-point energy
result = static_job(atoms)

# Geometry optimization
result = relax_job(atoms, relax_cell=True)

# Production relaxation (coarse then fine)
result = double_relax_flow(atoms)

# Access results
print(result.energy)
print(result.atoms)  # Relaxed structure
```

## Runners for Different Environments

```python
from vasp.runners import LocalRunner, MockRunner

# Local execution
runner = LocalRunner(command="mpirun -np 4 vasp_std")

# Testing without VASP
from vasp.runners import MockRunner, MockResults
mock = MockResults(energy=-10.5)
runner = MockRunner(results=mock)

# Use with calculator
calc = Vasp(atoms=atoms, runner=runner, ...)
```

## Key VASP Parameters

| Category | Parameters |
|----------|-----------|
| Basic | xc, encut, kpts, ismear, sigma |
| Relaxation | ibrion, isif, nsw, ediffg |
| Electronic | nelm, ediff, algo, lreal |
| Magnetic | ispin, magmom, lorbit |
| Output | lwave, lcharg, lvtot |

## Common Calculation Types

1. **Static**: Single-point energy
2. **Relax**: Geometry optimization (IBRION=2, ISIF=2 or 3)
3. **DOS**: Density of states (ISMEAR=-5, LORBIT=11)
4. **Bands**: Band structure (ICHARG=11, k-path)
5. **MD**: Molecular dynamics (IBRION=0, MDALGO)
6. **NEB**: Transition states (IMAGES, SPRING, LCLIMB)

## Documentation

For more details, see the vasp-ase documentation or run:
```bash
vasp-claude docs
```
"""

VASP_HELP_COMMAND = """Get help with VASP parameters and the vasp-ase interface.

Arguments: $ARGUMENTS

Provide detailed help based on the query:

1. **Parameter help** (e.g., "encut", "ispin"):
   - Explain what the parameter does
   - Show valid values and defaults
   - Give usage examples with vasp-ase

2. **Category help** (e.g., "relaxation", "magnetic"):
   - List all relevant parameters
   - Show typical settings
   - Provide code examples

3. **Workflow help** (e.g., "dos", "bands", "neb"):
   - Show complete workflow
   - Explain required parameters
   - Provide working code

4. **Troubleshooting** (e.g., "convergence", "memory"):
   - Common issues and solutions
   - Parameter adjustments

## Parameter Categories

- **electronic**: ENCUT, ISMEAR, SIGMA, EDIFF, NELM, ALGO
- **relaxation**: IBRION, ISIF, NSW, EDIFFG, POTIM
- **magnetic**: ISPIN, MAGMOM, LORBIT, LNONCOLLINEAR
- **hybrid**: LHFCALC, HFSCREEN, AEXX, ALGO
- **vdw**: IVDW, LVDW (use get_vdw_params())
- **dft+u**: LDAU, LDAUTYPE, LDAUL, LDAUU (use get_ldau_params())
- **md**: IBRION=0, MDALGO, TEBEG, TEEND, SMASS
- **neb**: IMAGES, SPRING, LCLIMB, ICHAIN

## Using Presets

```python
from vasp.parameters import (
    get_vdw_params,      # 'd2', 'd3', 'd3bj', 'ts'
    get_ldau_params,     # DFT+U setup
    get_hybrid_params,   # 'hse06', 'pbe0', 'b3lyp'
    get_soc_params,      # Spin-orbit coupling
    get_md_params,       # 'nvt', 'npt'
)
```
"""


def claude_install(args=None):
    """Install vasp-ase Claude Code skills globally."""
    claude_home = get_claude_home()
    commands_dir = claude_home / "commands"
    skills_dir = claude_home / "skills"

    # Create directories
    commands_dir.mkdir(parents=True, exist_ok=True)
    skills_dir.mkdir(parents=True, exist_ok=True)

    print(f"Installing vasp-ase Claude Code skills to {claude_home}")

    # Install global skill
    skill_file = skills_dir / "vasp.md"
    skill_file.write_text(GLOBAL_SKILL_CONTENT)
    print(f"  ✓ Installed skill: {skill_file}")

    # Install vasp-help command
    help_file = commands_dir / "vasp-help.md"
    help_file.write_text(VASP_HELP_COMMAND)
    print("  ✓ Installed command: /vasp-help")

    # Try to copy additional commands from package
    package_claude = get_package_claude_dir()
    if package_claude:
        package_commands = package_claude / "commands"
        if package_commands.exists():
            for cmd_file in package_commands.glob("*.md"):
                # Prefix with vasp- to avoid conflicts
                if not cmd_file.name.startswith("vasp-"):
                    target_name = f"vasp-{cmd_file.name}"
                else:
                    target_name = cmd_file.name
                target = commands_dir / target_name
                shutil.copy(cmd_file, target)
                print(f"  ✓ Installed command: /{target_name[:-3]}")

    print()
    print("Installation complete!")
    print()
    print("Available commands:")
    print("  /vasp-help <topic>  - Get VASP parameter help")
    print()
    print("The 'vasp' skill is now available globally.")
    print("Claude will automatically use it when you ask about VASP or DFT.")


def claude_uninstall(args=None):
    """Remove vasp-ase Claude Code skills."""
    claude_home = get_claude_home()

    removed = []

    # Remove skill
    skill_file = claude_home / "skills" / "vasp.md"
    if skill_file.exists():
        skill_file.unlink()
        removed.append(str(skill_file))

    # Remove commands with vasp- prefix
    commands_dir = claude_home / "commands"
    if commands_dir.exists():
        for cmd_file in commands_dir.glob("vasp-*.md"):
            cmd_file.unlink()
            removed.append(str(cmd_file))

    if removed:
        print("Removed vasp-ase Claude Code skills:")
        for f in removed:
            print(f"  - {f}")
    else:
        print("No vasp-ase Claude Code skills found to remove.")


def claude_status(args=None):
    """Show Claude Code skill installation status."""
    claude_home = get_claude_home()

    print(f"Claude Code home: {claude_home}")
    print()

    # Check skill
    skill_file = claude_home / "skills" / "vasp.md"
    if skill_file.exists():
        print("✓ vasp skill installed")
    else:
        print("✗ vasp skill not installed")

    # Check commands
    commands_dir = claude_home / "commands"
    vasp_commands = list(commands_dir.glob("vasp-*.md")) if commands_dir.exists() else []

    if vasp_commands:
        print(f"✓ {len(vasp_commands)} vasp commands installed:")
        for cmd in vasp_commands:
            print(f"    /{cmd.stem}")
    else:
        print("✗ No vasp commands installed")

    print()
    print("Run 'vasp-claude install' to install or update skills.")


def claude_docs(args=None):
    """Show documentation location and open if possible."""
    try:
        import vasp

        package_dir = Path(vasp.__file__).parent.parent
        docs_dir = package_dir / "docs"

        if docs_dir.exists():
            print(f"Documentation source: {docs_dir}")

            # Check for built docs
            built_docs = docs_dir / "_build" / "html" / "index.html"
            if built_docs.exists():
                print(f"Built documentation: {built_docs}")
                print()
                print("To open in browser:")
                print(f"  open {built_docs}  # macOS")
                print(f"  xdg-open {built_docs}  # Linux")
            else:
                print()
                print("Documentation not built. To build:")
                print("  pip install jupyter-book")
                print(f"  jupyter-book build {docs_dir}")
        else:
            print("Documentation not found in package.")
            print("Visit: https://github.com/jkitchin/vasp")

    except ImportError:
        print("vasp-ase package not found.")
        print("Install with: pip install vasp-ase")


def status_command():
    """Show VASP calculation status (placeholder for job monitoring)."""
    print("VASP-ASE Status")
    print("=" * 40)

    try:
        import vasp

        print(f"Package version: {vasp.__version__}")
    except ImportError:
        print("Package not installed")
        return

    # Check environment
    pp_path = os.environ.get("VASP_PP_PATH", "Not set")
    vasp_cmd = os.environ.get("ASE_VASP_COMMAND", "Not set")

    print(f"VASP_PP_PATH: {pp_path}")
    print(f"ASE_VASP_COMMAND: {vasp_cmd}")


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="vasp-claude", description="VASP-ASE command-line tools and Claude Code integration"
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Install command
    install_parser = subparsers.add_parser("install", help="Install Claude Code skills globally")
    install_parser.set_defaults(func=claude_install)

    # Uninstall command
    uninstall_parser = subparsers.add_parser("uninstall", help="Remove Claude Code skills")
    uninstall_parser.set_defaults(func=claude_uninstall)

    # Status command
    status_parser = subparsers.add_parser("status", help="Show installation status")
    status_parser.set_defaults(func=claude_status)

    # Docs command
    docs_parser = subparsers.add_parser("docs", help="Show documentation location")
    docs_parser.set_defaults(func=claude_docs)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        print()
        print("Quick start:")
        print("  vasp-claude install   # Install Claude Code skills")
        print("  vasp-claude status    # Check installation")
        print("  vasp-claude docs      # Find documentation")
        sys.exit(0)

    args.func(args)


def vaspsum():
    """Summarize VASP calculations.

    CLI tool to display quick summaries of VASP calculation results.
    """
    import argparse

    from vasp import Vasp
    from vasp.exceptions import VaspException

    def format_energy(energy, natoms):
        if energy is None:
            return "Energy: Not available"
        return f"Energy: {energy:.6f} eV ({energy/natoms:.6f} eV/atom)"

    def format_forces(forces):
        if forces is None:
            return "Forces: Not available"
        fmax = (forces**2).sum(axis=1).max() ** 0.5
        return f"Max force: {fmax:.6f} eV/Å"

    def format_stress(stress):
        if stress is None:
            return "Stress: Not available"
        stress_gpa = stress * 0.1
        pressure = -stress_gpa[:3].mean()
        return f"Pressure: {pressure:.3f} GPa"

    def print_summary(calc, verbose=False):
        try:
            atoms = calc.atoms
            if atoms is None:
                # Try to read from CONTCAR/POSCAR
                try:
                    atoms = calc.load_atoms()
                    calc.atoms = atoms
                except Exception:
                    pass
            if atoms is None:
                print("  No structure found")
                return

            print(f"  Formula: {atoms.get_chemical_formula()}")
            print(f"  Atoms: {len(atoms)}")

            try:
                energy = calc.results.get("energy")
                forces = calc.results.get("forces")
                stress = calc.results.get("stress")

                print(f"  {format_energy(energy, len(atoms))}")
                if verbose and forces is not None:
                    print(f"  {format_forces(forces)}")
                if verbose and stress is not None:
                    print(f"  {format_stress(stress)}")

                converged = calc.results.get("converged", None)
                if converged is not None:
                    status = "Converged" if converged else "NOT CONVERGED"
                    print(f"  Status: {status}")

            except (VaspException, FileNotFoundError) as e:
                print(f"  Results not available: {e}")

        except Exception as e:
            print(f"  Error reading calculation: {e}")

    def print_parameters(calc, verbose=False):
        print("\nParameters:")
        params = calc.parameters
        key_params = ["xc", "encut", "kpts", "ismear", "sigma"]

        for key in key_params:
            if key in params:
                print(f"  {key}: {params[key]}")

        if verbose:
            print("\nAll parameters:")
            for key, value in sorted(params.items()):
                if key not in key_params:
                    print(f"  {key}: {value}")

    def print_json(calc, pretty=False):
        data = {
            "directory": calc.directory,
            "parameters": calc.parameters,
            "results": calc.results,
        }

        if calc.atoms:
            data["formula"] = calc.atoms.get_chemical_formula()
            data["natoms"] = len(calc.atoms)
            data["positions"] = calc.atoms.positions.tolist()
            data["cell"] = calc.atoms.cell.tolist()
            data["symbols"] = calc.atoms.get_chemical_symbols()

        if pretty:
            print(json.dumps(data, indent=2, default=str))
        else:
            print(json.dumps(data, default=str))

    parser = argparse.ArgumentParser(
        description="Summarize VASP calculations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  vaspsum                    # Summarize current directory
  vaspsum calc1 calc2        # Summarize multiple directories
  vaspsum -v .               # Verbose output
  vaspsum --json calc1       # JSON output
  vaspsum --view .           # View structure
        """,
    )

    parser.add_argument(
        "dirs",
        nargs="*",
        default=["."],
        help="Directories to summarize (default: current directory)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Show detailed information")
    parser.add_argument(
        "--params", "--describe", action="store_true", help="Show calculation parameters"
    )
    parser.add_argument("--json", action="store_true", help="Output in JSON format")
    parser.add_argument(
        "--json-pretty",
        "--jsonpp",
        action="store_true",
        help="Output in pretty-printed JSON format",
    )
    parser.add_argument("--view", "-p", action="store_true", help="View final structure")

    args = parser.parse_args()

    for directory in args.dirs:
        if not os.path.isdir(directory):
            print(f"Error: {directory} does not exist!", file=sys.stderr)
            continue

        print(f"\n{directory}")
        print("=" * 60)

        try:
            calc = Vasp(directory)

            # Try to load atoms and results from existing files
            try:
                calc.atoms = calc.load_atoms()
                calc._setup_sorting(calc.atoms)
            except Exception:
                pass
            try:
                calc.read_results()
            except Exception:
                pass

            if args.json or args.json_pretty:
                print_json(calc, pretty=args.json_pretty)
            elif args.params:
                print_summary(calc, verbose=args.verbose)
                print_parameters(calc, verbose=args.verbose)
            else:
                print_summary(calc, verbose=args.verbose)

            if args.view and calc.atoms:
                from ase.visualize import view as ase_view

                ase_view(calc.atoms)

        except Exception as e:
            print(f"Error processing {directory}: {e}", file=sys.stderr)
            if args.verbose:
                import traceback

                traceback.print_exc()


def vasp_debug():
    """Output comprehensive diagnostic information for VASP setup.

    This tool helps troubleshoot VASP installation and configuration issues
    by displaying system info, environment variables, library versions,
    executable paths, and pseudopotential availability.
    """
    import platform
    import subprocess

    def section(title):
        print()
        print("=" * 60)
        print(f" {title}")
        print("=" * 60)

    def check_command(cmd, description=None):
        """Check if a command exists and return its path."""
        path = shutil.which(cmd)
        if path:
            return f"✓ {cmd}: {path}"
        return f"✗ {cmd}: not found"

    def get_version(cmd):
        """Try to get version from a command."""
        try:
            result = subprocess.run([cmd, "--version"], capture_output=True, text=True, timeout=5)
            return result.stdout.strip().split("\n")[0]
        except Exception:
            return None

    def check_package(name, import_name=None):
        """Check if a Python package is installed and get version."""
        import_name = import_name or name
        try:
            mod = __import__(import_name)
            version = getattr(mod, "__version__", "unknown")
            return f"✓ {name}: {version}"
        except ImportError:
            return f"✗ {name}: not installed"

    def check_pp_directory(base_path, functional="potpaw_PBE"):
        """Check pseudopotential directory structure."""
        pp_path = Path(base_path) / functional
        if not pp_path.exists():
            return None, []

        elements = []
        for item in sorted(pp_path.iterdir()):
            if item.is_dir():
                potcar = item / "POTCAR"
                if potcar.exists():
                    elements.append(item.name)

        return str(pp_path), elements

    # Header
    print()
    print("VASP-ASE Debug Information")
    print("Generated:", __import__("datetime").datetime.now().isoformat())

    # System Information
    section("System Information")
    print(f"Platform:        {platform.platform()}")
    print(f"Architecture:    {platform.machine()}")
    print(f"Python:          {platform.python_version()}")
    print(f"Python path:     {sys.executable}")
    print(f"Working dir:     {os.getcwd()}")

    # Environment Variables
    section("Environment Variables (VASP-related)")

    vasp_env_vars = [
        "VASP_PP_PATH",
        "ASE_VASP_COMMAND",
        "VASP_COMMAND",
        "VASP_SCRIPT",
        "ASE_VASP_VDW",
        "VASP_PP",
        "POTCAR_PATH",
    ]

    for var in vasp_env_vars:
        value = os.environ.get(var)
        if value:
            print(f"✓ {var}={value}")
        else:
            print(f"  {var}: not set")

    # Additional useful environment variables
    print()
    print("Other relevant environment variables:")
    other_vars = ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "PATH", "LD_LIBRARY_PATH"]
    for var in other_vars:
        value = os.environ.get(var)
        if value:
            if var in ["PATH", "LD_LIBRARY_PATH"]:
                # Show first few paths
                paths = value.split(":")[:5]
                print(f"  {var}:")
                for p in paths:
                    print(f"    {p}")
                if len(value.split(":")) > 5:
                    print(f"    ... and {len(value.split(':')) - 5} more")
            else:
                print(f"  {var}={value}")

    # Python Packages
    section("Python Packages")

    packages = [
        ("vasp-ase", "vasp"),
        ("ase", "ase"),
        ("numpy", "numpy"),
        ("scipy", "scipy"),
        ("matplotlib", "matplotlib"),
        ("spglib", "spglib"),
        ("phonopy", "phonopy"),
        ("pymatgen", "pymatgen"),
        ("mp-api", "mp_api"),
        ("icet", "icet"),
        ("mpi4py", "mpi4py"),
    ]

    for name, import_name in packages:
        print(check_package(name, import_name))

    # VASP Executables
    section("VASP Executables")

    vasp_commands = [
        "vasp_std",
        "vasp_gam",
        "vasp_ncl",
        "vasp",
    ]

    for cmd in vasp_commands:
        print(check_command(cmd))

    # Check ASE_VASP_COMMAND specifically
    ase_cmd = os.environ.get("ASE_VASP_COMMAND")
    if ase_cmd:
        # Extract the actual vasp command (might have mpirun prefix)
        parts = ase_cmd.split()
        vasp_part = parts[-1] if parts else None
        if vasp_part:
            path = shutil.which(vasp_part)
            if path:
                print(f"✓ ASE_VASP_COMMAND resolves to: {path}")
            else:
                print(f"✗ ASE_VASP_COMMAND binary not found: {vasp_part}")

    # MPI Configuration
    section("MPI Configuration")

    mpi_commands = ["mpirun", "mpiexec", "srun", "aprun", "orterun"]
    for cmd in mpi_commands:
        result = check_command(cmd)
        if "✓" in result:
            print(result)
            version = get_version(cmd)
            if version:
                print(f"    Version: {version}")

    # Check for common MPI implementations
    mpi_libs = ["openmpi", "mpich", "intelmpi"]
    print()
    print("MPI implementation hints:")
    for lib in mpi_libs:
        # Check if library directory exists or module loaded
        for path_dir in os.environ.get("PATH", "").split(":"):
            if lib in path_dir.lower():
                print(f"  Possible {lib}: {path_dir}")
                break

    # Pseudopotentials
    section("Pseudopotentials")

    pp_path = os.environ.get("VASP_PP_PATH")
    if not pp_path:
        print("✗ VASP_PP_PATH not set")
        print()
        print("To set up pseudopotentials:")
        print("  export VASP_PP_PATH=/path/to/vasp/potpaw")
        print()
        print("Expected directory structure:")
        print("  $VASP_PP_PATH/")
        print("  ├── potpaw_PBE/")
        print("  │   ├── Ag/POTCAR")
        print("  │   ├── Al/POTCAR")
        print("  │   └── ...")
        print("  ├── potpaw_LDA/")
        print("  └── potpaw_GGA/")
    else:
        pp_base = Path(pp_path)
        if not pp_base.exists():
            print(f"✗ VASP_PP_PATH directory does not exist: {pp_path}")
        else:
            print(f"✓ VASP_PP_PATH: {pp_path}")
            print()

            # Check for different functionals
            functionals = ["potpaw_PBE", "potpaw_LDA", "potpaw_GGA", "PBE", "LDA", "GGA"]
            found_any = False

            for func in functionals:
                func_path, elements = check_pp_directory(pp_path, func)
                if func_path:
                    found_any = True
                    print(f"✓ {func}: {len(elements)} elements")

                    # Show some sample elements
                    if elements:
                        sample = elements[:10]
                        print(f"    Sample: {', '.join(sample)}")
                        if len(elements) > 10:
                            print(f"    ... and {len(elements) - 10} more")

                        # Check for common variants
                        variants = []
                        for elem in elements[:20]:
                            if "_" in elem:
                                base, variant = elem.rsplit("_", 1)
                                if variant not in variants:
                                    variants.append(variant)
                        if variants:
                            print(f"    Variants found: {', '.join(variants[:5])}")

            if not found_any:
                print("✗ No pseudopotential directories found")
                print("  Expected: potpaw_PBE, potpaw_LDA, etc.")

    # ASE Configuration
    section("ASE Configuration")

    try:
        import ase

        print("✓ ASE imported successfully")

        # Check ASE data directory
        ase_dir = Path(ase.__file__).parent
        print(f"  ASE location: {ase_dir}")

        # Check for ASE database
        db_path = Path.home() / ".ase"
        if db_path.exists():
            print(f"  ASE config dir: {db_path}")

    except ImportError as e:
        print(f"✗ ASE import error: {e}")

    # vasp-ase specific
    section("vasp-ase Package")

    try:
        import vasp

        print(f"✓ Version: {vasp.__version__}")
        print(f"  Location: {Path(vasp.__file__).parent}")

        # Check for runners
        from vasp import runners

        available_runners = []
        for name in ["LocalRunner", "MockRunner", "SlurmRunner", "InteractiveRunner"]:
            if hasattr(runners, name):
                available_runners.append(name)
        print(f"  Runners: {', '.join(available_runners)}")

        # Check for recipes
        try:
            from vasp import recipes  # noqa: F401

            print("  Recipes: available")
        except ImportError:
            print("  Recipes: not available")

    except ImportError as e:
        print(f"✗ vasp-ase import error: {e}")

    # Quick connectivity test
    section("Quick Tests")

    print("Calculator instantiation test:")
    try:
        import numpy as np
        from ase.build import bulk

        from vasp import Vasp
        from vasp.runners import MockResults, MockRunner

        atoms = bulk("Si", "diamond", a=5.43)
        mock = MockResults(energy=-10.0, forces=np.zeros((2, 3)))
        calc = Vasp(
            atoms=atoms,
            runner=MockRunner(results=mock),
            xc="PBE",
            encut=300,
        )
        energy = calc.potential_energy
        print(f"✓ Mock calculation successful: E = {energy:.4f} eV")
    except Exception as e:
        print(f"✗ Calculator test failed: {e}")

    print()
    print("POTCAR path test:")
    pp_path = os.environ.get("VASP_PP_PATH")
    if pp_path:
        # Check if Si POTCAR exists
        si_potcar = Path(pp_path) / "potpaw_PBE" / "Si" / "POTCAR"
        if si_potcar.exists():
            print(f"✓ Si POTCAR found: {si_potcar}")
        else:
            # Try alternate path
            si_potcar = Path(pp_path) / "PBE" / "Si" / "POTCAR"
            if si_potcar.exists():
                print(f"✓ Si POTCAR found: {si_potcar}")
            else:
                print("✗ Si POTCAR not found in expected locations")
    else:
        print("✗ Cannot test POTCAR - VASP_PP_PATH not set")

    # Summary
    section("Summary")

    issues = []

    if not os.environ.get("VASP_PP_PATH"):
        issues.append("VASP_PP_PATH not set - cannot generate POTCAR files")

    if not os.environ.get("ASE_VASP_COMMAND"):
        issues.append("ASE_VASP_COMMAND not set - cannot run VASP automatically")

    vasp_found = any(shutil.which(cmd) for cmd in vasp_commands)
    if not vasp_found:
        issues.append("No VASP executable found in PATH")

    if issues:
        print("Issues found:")
        for issue in issues:
            print(f"  ⚠ {issue}")
    else:
        print("✓ No critical issues detected")

    print()
    print("For more help, see: https://kitchingroup.cheme.cmu.edu/vasp/")
    print()


if __name__ == "__main__":
    main()
