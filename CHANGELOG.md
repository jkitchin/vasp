# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Makefile for development workflow (test, lint, build, publish)
- CONTRIBUTING.md with development and publishing guidelines
- CODE_REVIEW.md with comprehensive code analysis
- MIGRATION.md guide for v1.x to v2.0 migration

### Changed
- Fixed `.claude/settings.json` to use new hooks format
- Updated `bin/vaspsum` to use modern vasp-ase API
- Updated `vaspy-mode.el` with hardcoded VASP parameters (removed dependency on deleted modules)
- Improved error messages with helpful suggestions and documentation links

### Removed
- Legacy .pyc files (22 files from 2016-2017)
- `dodo.py` (deprecated build tool)
- `bin/runvasp.py` (depended on deleted vasprc module)
- `__pycache__` directories

## [2.0.0] - 2025-12-27

**Major rewrite with breaking changes.** This release represents a complete modernization of the package.

### Added

#### Core Features
- Modern ASE calculator interface with full type hints
- Parameter preset functions (`get_vdw_params()`, `get_ldau_params()`, `get_hybrid_params()`)
- Pluggable execution backends (LocalRunner, SlurmRunner, KubernetesRunner)
- InteractiveRunner for persistent VASP sessions with socket I/O
- Socket I/O implementation (i-PI protocol) for remote VASP communication
- Exception-based async workflow pattern (VaspSubmitted, VaspQueued, VaspRunning)
- Workflow recipes with @job and @flow decorators (quacc-style)
- Vector database with per-atom embeddings using libSQL
- Comprehensive mixin architecture (IOMixin, ElectronicMixin, AnalysisMixin, DynamicsMixin)

#### Documentation and Examples
- 21 progressive tutorials (beginner → intermediate → advanced)
  - 01: Getting Started
  - 02: Convergence Testing
  - 03: Relaxation
  - 04: Equation of State
  - 05: Density of States
  - 06: Band Structure
  - 07: Magnetism
  - 08: Surfaces
  - 09: Adsorption
  - 10: Reactions
  - 11: Phonons
  - 12: DFT+U
  - 13: Hybrid Functionals
  - 14: Van der Waals
  - 15: Workflows
  - 16: NEB
  - 17: Vibrations
  - 18: 3D Visualization
  - 19: Pseudopotentials
  - 20: Interactive Mode
  - 21: Cluster Expansion (with ICET)
- Jupyter Book documentation structure
- CLAUDE.md for AI assistant integration
- Comprehensive API documentation in docstrings

#### Claude Code Integration
- 13 project-local slash commands:
  - `/docs` - Open documentation
  - `/examples` - List all examples
  - `/tutorial <n>` - View specific tutorial
  - `/test` - Run test suite
  - `/new-example` - Create new example
  - `/architecture` - Review codebase
  - `/lint` - Run code quality checks
  - `/build-docs` - Build Jupyter Book
  - `/status` - Project status
  - `/fix-job <dir>` - Diagnose and fix failed job
  - `/watch-job <dir>` - Monitor VASP job
  - `/vasp-help <topic>` - Parameter reference
  - `/notebooks` - Explore notebooks
- 3 specialized skills:
  - `vasp-calculation` - Creating VASP calculations
  - `job-watcher` - Monitoring jobs
  - `troubleshoot` - Error diagnosis
- Global skill installer via `vasp-claude` CLI
- SessionStart hook showing package version

#### Testing and Quality
- 282 comprehensive tests
- MockRunner for testing without VASP installation
- Test fixtures and shared utilities
- Coverage reporting integration
- Type hints throughout entire codebase
- Ruff linting configuration
- mypy type checking

#### Command-Line Tools
- `vasp-claude install` - Install global Claude Code skills
- `vasp-claude status` - Check installation status
- `vasp-claude uninstall` - Remove global skills
- `vasp-claude docs` - Open documentation
- `vasp-status` - Show VASP environment status
- `vaspsum` - Summarize VASP calculations

#### Parameter Support
- Van der Waals corrections: D2, D3, D3-BJ, TS, MBD, VdW-DF variants
- DFT+U with HubbardU dataclass and common element defaults
- Hybrid functionals: HSE06, PBE0, B3LYP, HSE03
- Spin-orbit coupling parameters
- Machine learning force fields (MLFF)
- Molecular dynamics presets (NVE, NVT, NPT)
- Phonon calculation presets

#### Recipes and Workflows
- `static_job()` - Single-point energy calculation
- `relax_job()` - Geometry optimization
- `double_relax_flow()` - Two-stage relaxation (coarse → fine)
- `slab_static_job()` - Surface static calculation
- `slab_relax_job()` - Surface relaxation
- `phonon_flow()` - Phonon calculation workflow
- Recipe decorators: `@job`, `@flow`, `@subflow`
- Workflow engine integration (Prefect, Dask, Parsl, Covalent)

#### Analysis Tools
- Bader charge analysis
- Elastic moduli extraction
- Band structure calculation
- Density of states
- Charge density analysis
- Molecular dynamics trajectory analysis

### Changed
- Complete rewrite of calculator class using ASE's modern Calculator interface
- Switched from setup.py to pyproject.toml (PEP 621)
- Changed package structure to use mixins instead of monkeypatching
- Updated all examples to use new API
- Modernized test suite with pytest
- Replaced nosetests/pep8/pyflakes with pytest/ruff/mypy
- Changed from README.org to README.md
- Updated license to MIT

### Removed

#### Legacy Code (Breaking Changes)
- **Monkeypatch-based calculator implementation**
- **Old Vasp class with get/set methods**
- MongoDB integration (`vasp.mongo` module)
- Legacy validation module (`vasp.validate`)
- Old reader/writer modules
- Legacy runner system
- Old test suite based on nosetests
- `vasp.vasprc` configuration module
- All getters/setters modules
- `vasp.monkeypatch` module
- `vasp.serialize` module
- Legacy exception classes
- Old bader/bandstructure/elastic_moduli modules (functionality moved to mixins)

#### Build and Configuration
- setup.py (replaced by pyproject.toml)
- requirements.txt (replaced by pyproject.toml dependencies)
- README.org (replaced by README.md)
- .travis.yml (likely replaced by GitHub Actions)

### Migration Notes

**This is a major breaking release.** Code written for v1.x will need updates. See MIGRATION.md for detailed migration instructions.

Key changes requiring code updates:
1. Calculator initialization now requires `atoms` parameter
2. `calc.calculate()` replaced by accessing properties (e.g., `calc.potential_energy`)
3. Parameters passed as keyword arguments instead of via `set()` method
4. No more MongoDB integration
5. Different import paths for some functionality

### Fixed
- Numerous bug fixes and stability improvements
- Better error handling and exceptions
- Improved POTCAR file handling
- Fixed stress tensor handling in socket I/O
- Better convergence detection

### Security
- Removed hardcoded paths
- Better environment variable handling
- Improved input validation

## [1.x.x] - Legacy

Previous versions used a monkeypatch-based approach. See git history for details on v1.x releases.

### Notable v1.x Features (Now Removed)
- Monkeypatch-based ASE calculator modification
- MongoDB integration for results storage
- Dynamic keyword validation
- Object-oriented getter/setter interface
- README.org with org-mode examples

---

## Version History Summary

- **2.0.0** (2025-12-27): Complete rewrite with modern architecture
- **1.x.x**: Legacy monkeypatch-based implementation

---

## How to Upgrade

See [MIGRATION.md](MIGRATION.md) for detailed upgrade instructions from v1.x to v2.0.

## Links

- [Repository](https://github.com/jkitchin/vasp)
- [Documentation](https://github.com/jkitchin/vasp/tree/master/docs)
- [Issue Tracker](https://github.com/jkitchin/vasp/issues)
- [PyPI](https://pypi.org/project/vasp-ase/)
