Run code quality checks on the project.

Execute the following linting and type checking:

1. **Ruff (linting)**
   ```bash
   ruff check vasp/
   ```

2. **Ruff (formatting check)**
   ```bash
   ruff format --check vasp/
   ```

3. **MyPy (type checking)**
   ```bash
   mypy vasp/ --ignore-missing-imports
   ```

Report any issues found and offer to fix them:
- For ruff: `ruff check --fix vasp/`
- For formatting: `ruff format vasp/`

Summarize the code quality status.
