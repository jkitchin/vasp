Build and preview the Jupyter Book documentation.

Execute the documentation build process:

1. Check if jupyter-book is installed:
   ```bash
   pip show jupyter-book
   ```

2. If not installed, install documentation dependencies:
   ```bash
   pip install -r docs/requirements.txt
   ```

3. Build the documentation:
   ```bash
   jupyter-book build docs
   ```

4. Report build status and any warnings/errors

5. Show how to view the built documentation:
   - Local path: docs/_build/html/index.html
   - Suggest opening in browser

6. If there are build errors:
   - Analyze the error messages
   - Suggest fixes
   - Offer to fix automatically if appropriate
