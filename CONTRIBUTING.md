# Contributing to VASP-ASE

Thank you for your interest in contributing to VASP-ASE!

## Development Setup

1. **Clone the repository:**
   ```bash
   git clone https://github.com/jkitchin/vasp.git
   cd vasp
   ```

2. **Install in development mode:**
   ```bash
   make install-dev
   # or manually:
   pip install -e ".[dev]"
   ```

## Development Workflow

The project uses a `Makefile` for common development tasks:

### Running Tests

```bash
# Run all tests
make test

# Run tests with coverage report
make test-coverage
```

### Code Quality

```bash
# Run linter (ruff)
make lint

# Run type checker (mypy)
make typecheck

# Auto-format code
make format

# Run all checks (lint + typecheck + test)
make check
```

### Building and Publishing

```bash
# Build package
make build

# Publish to TestPyPI (for testing)
make publish-test

# Publish to PyPI (production)
make publish
```

### Documentation

```bash
# Build Jupyter Book documentation
make docs-build

# Build and serve locally at http://localhost:8000
make docs-serve
```

### Cleaning

```bash
# Remove build artifacts
make clean

# Remove all generated files including cache
make clean-all
```

### Quick Status

```bash
# Show version, git status, and uncommitted changes
make status
```

## Publishing a New Release

1. **Update the version** in `pyproject.toml`:
   ```toml
   version = "2.0.1"  # Update this
   ```

2. **Update CHANGELOG** (if applicable):
   - Add release notes for the new version

3. **Commit changes:**
   ```bash
   git add pyproject.toml CHANGELOG.md
   git commit -m "Bump version to 2.0.1"
   git tag v2.0.1
   git push origin master --tags
   ```

4. **Test the build locally:**
   ```bash
   make check  # Run all tests and checks
   ```

5. **Test on TestPyPI first** (optional but recommended):
   ```bash
   make publish-test
   # Then test installation:
   pip install --index-url https://test.pypi.org/simple/ vasp-ase
   ```

6. **Publish to PyPI:**
   ```bash
   make publish
   ```

   You'll need PyPI credentials. Configure them with:
   ```bash
   # Create ~/.pypirc with:
   [distutils]
   index-servers =
       pypi
       testpypi

   [pypi]
   username = __token__
   password = pypi-...  # Your PyPI API token

   [testpypi]
   repository = https://test.pypi.org/legacy/
   username = __token__
   password = pypi-...  # Your TestPyPI API token
   ```

   Or use environment variables:
   ```bash
   export TWINE_USERNAME=__token__
   export TWINE_PASSWORD=pypi-...
   ```

## Code Style

- Follow PEP 8 style guide
- Use type hints for all function signatures
- Maximum line length: 100 characters
- Run `make format` before committing

## Testing

- Write tests for all new features
- Maintain test coverage above 80%
- Use `MockRunner` for testing without VASP installed
- Place tests in `vasp/tests/`

## Documentation

- Update docstrings for any modified functions
- Add examples to docstrings when helpful
- Update Jupyter Book tutorials for new features
- Follow Google-style docstrings

## Commit Messages

Use clear, descriptive commit messages:

```
Add support for MLFF calculations

- Implement ML_LMLFF parameter handling
- Add test cases for MLFF workflows
- Update documentation with MLFF example
```

## Questions?

Open an issue on GitHub or contact the maintainers.
