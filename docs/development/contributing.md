# Contributing

We welcome contributions to the vasp-ase project!

## Development Setup

```bash
git clone https://github.com/jkitchin/vasp.git
cd vasp
pip install -e ".[dev]"
```

## Running Tests

```bash
pytest
```

With coverage:

```bash
pytest --cov=vasp
```

## Code Style

We use:
- **ruff** for linting
- **mypy** for type checking

```bash
ruff check .
mypy vasp
```

## Pull Request Process

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Run tests and linting
6. Submit a pull request

## Documentation

Build docs locally:

```bash
pip install -r docs/requirements.txt
jupyter-book build docs
```

## Adding Examples

Examples go in `examples/XX_topic/`:

```
examples/
├── XX_topic/
│   ├── run.py       # Main script
│   └── README.md    # Explanation
```

## Code of Conduct

Be respectful and inclusive in all interactions.
