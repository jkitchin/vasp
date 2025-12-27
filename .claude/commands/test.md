Run the project test suite.

Execute the following:

1. Run pytest with verbose output:
   ```bash
   pytest -v
   ```

2. If tests fail, analyze the failures and suggest fixes.

3. Optionally run with coverage:
   ```bash
   pytest --cov=vasp --cov-report=term-missing
   ```

Report:
- Total tests passed/failed/skipped
- Any test failures with details
- Coverage summary if requested
