Show the current project status including tests, git, and package info.

Run the following checks and provide a summary:

1. **Package Status**
   ```bash
   python -c "import vasp; print(f'Version: {vasp.__version__}')"
   ```

2. **Git Status**
   ```bash
   git status --short
   git log -1 --oneline
   ```

3. **Test Status** (quick check)
   ```bash
   pytest --collect-only -q 2>/dev/null | tail -5
   ```

4. **Dependencies**
   ```bash
   pip show ase numpy scipy | grep -E "^(Name|Version):"
   ```

5. **Environment**
   - Check if VASP_PP_PATH is set
   - Check if ASE_VASP_COMMAND is set

Provide a concise summary of the project state.
