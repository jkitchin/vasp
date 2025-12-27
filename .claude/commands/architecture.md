Review the current codebase structure and architecture.

Perform a comprehensive review:

1. **Package Structure**
   - Read vasp/__init__.py to understand exports
   - List all modules in vasp/
   - Identify key classes and functions

2. **Calculator Architecture**
   - Show the Vasp class interface
   - Explain the runner abstraction
   - Document available runners

3. **Recipe System**
   - List all available recipes in vasp/recipes/
   - Explain the @job and @flow decorators
   - Show workflow engine integration

4. **Parameter Presets**
   - List all preset functions in vasp/parameters.py
   - Show available VdW methods
   - Document DFT+U setup

5. **Test Coverage**
   - List test files
   - Show test categories
   - Identify any gaps

Provide a summary suitable for a developer wanting to contribute.
