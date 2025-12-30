.PHONY: help install install-dev test test-coverage lint typecheck format clean clean-all build publish publish-test docs docs-build docs-serve tutorial tutorials tutorials-clean

# Default target
help:
	@echo "VASP-ASE Development Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  install        Install package in editable mode"
	@echo "  install-dev    Install package with dev dependencies"
	@echo "  test           Run tests with pytest"
	@echo "  test-coverage  Run tests with coverage report"
	@echo "  lint           Run ruff linter"
	@echo "  typecheck      Run mypy type checker"
	@echo "  format         Format code with ruff"
	@echo "  check          Run all checks (lint, typecheck, test)"
	@echo "  clean          Remove build artifacts"
	@echo "  clean-all      Remove all generated files including cache"
	@echo "  build          Build package for distribution"
	@echo "  publish-test   Publish to TestPyPI"
	@echo "  publish        Publish to PyPI"
	@echo "  docs-build     Build Jupyter Book documentation"
	@echo "  docs-serve     Build and serve documentation locally"
	@echo "  tutorial T=x   Run a single tutorial (by name or number)"
	@echo "  tutorials      Run all 21 tutorials"
	@echo "  tutorials-clean Clear notebook outputs"
	@echo ""

# Installation targets
install:
	pip install -e .

install-dev:
	pip install -e ".[dev]"

# Testing targets
test:
	pytest -v

test-coverage:
	pytest --cov=vasp --cov-report=term-missing --cov-report=html
	@echo ""
	@echo "Coverage report generated in htmlcov/index.html"

# Code quality targets
lint:
	ruff check .

typecheck:
	mypy vasp

format:
	ruff check --fix .
	ruff format .

# Run all checks
check: lint typecheck test

# Cleaning targets
clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf .eggs/
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name '*.pyc' -delete
	find . -type f -name '*.pyo' -delete
	find . -type f -name '*~' -delete

clean-all: clean
	rm -rf .pytest_cache/
	rm -rf .coverage
	rm -rf htmlcov/
	rm -rf .mypy_cache/
	rm -rf .ruff_cache/
	rm -rf docs/_build/

# Build targets
build: clean
	python -m build

# Publishing targets
publish-test: build
	@echo "Publishing to TestPyPI..."
	@echo "You will need TestPyPI credentials."
	python -m twine upload --repository testpypi dist/*
	@echo ""
	@echo "Test installation with:"
	@echo "  pip install --index-url https://test.pypi.org/simple/ vasp-ase"

publish: build
	@echo "WARNING: This will publish to PyPI!"
	@echo "Press Ctrl+C to cancel, Enter to continue..."
	@read dummy
	python -m twine upload dist/*
	@echo ""
	@echo "Published to PyPI!"
	@echo "Install with: pip install vasp-ase"

# Documentation targets
docs-build:
	jupyter-book build docs

docs-serve: docs-build
	@echo "Serving documentation at http://localhost:8000"
	@cd docs/_build/html && python -m http.server 8000

# Version management
version:
	@grep '^version = ' pyproject.toml | sed 's/version = "\(.*\)"/\1/'

# Quick status check
status:
	@echo "VASP-ASE Development Status"
	@echo "============================"
	@echo "Version: $$(grep '^version = ' pyproject.toml | sed 's/version = "\(.*\)"/\1/')"
	@echo ""
	@echo "Git Status:"
	@git status --short
	@echo ""
	@echo "Git Branch:"
	@git branch --show-current
	@echo ""
	@echo "Uncommitted changes:"
	@git diff --shortstat

# Tutorial targets (notebooks are in docs/tutorials/)
NOTEBOOKS := $(sort $(wildcard docs/tutorials/*.ipynb))

tutorials:
	@echo "Running all tutorial notebooks..."
	@for nb in $(NOTEBOOKS); do \
		name=$$(basename $$nb .ipynb); \
		echo "  [$$name] Running..."; \
		jupyter nbconvert --to notebook --execute --inplace "$$nb" || exit 1; \
	done
	@echo "All tutorials completed."

# Run a single tutorial: make tutorial T=01 or T=getting_started
tutorial:
ifdef T
	@nb=$$(ls docs/tutorials/$(T)*.ipynb docs/tutorials/*$(T)*.ipynb 2>/dev/null | head -1); \
	if [ -z "$$nb" ]; then \
		echo "Error: No tutorial found matching '$(T)'"; \
		echo ""; \
		echo "Available tutorials:"; \
		ls docs/tutorials/*.ipynb 2>/dev/null | xargs -n1 basename | sed 's/.ipynb//'; \
		exit 1; \
	fi; \
	echo "Running $$nb..."; \
	jupyter nbconvert --to notebook --execute --inplace "$$nb"
else
	@echo "Usage: make tutorial T=<name>"
	@echo ""
	@echo "Examples:"
	@echo "  make tutorial T=01"
	@echo "  make tutorial T=getting_started"
	@echo "  make tutorial T=01_getting_started"
	@echo ""
	@echo "Available tutorials:"; \
	ls docs/tutorials/*.ipynb 2>/dev/null | xargs -n1 basename | sed 's/.ipynb//'
endif

tutorials-clean:
	@echo "Clearing notebook outputs..."
	@for nb in $(NOTEBOOKS); do \
		echo "  Clearing $$nb..."; \
		jupyter nbconvert --clear-output --inplace "$$nb"; \
	done
	@echo "Notebook outputs cleared."
