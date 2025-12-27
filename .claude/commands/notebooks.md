Open and explore Jupyter notebooks in this project.

Search for and list all Jupyter notebooks (.ipynb files) in the project:

1. Use glob to find all .ipynb files
2. For each notebook found:
   - Show its location
   - Read and summarize its contents
   - List the cell types (code/markdown)

3. If no notebooks exist, offer to:
   - Convert an example run.py to a notebook
   - Create a new tutorial notebook

4. If notebooks exist, ask which one the user wants to explore or run.

For running notebooks, suggest:
```bash
jupyter notebook <path>
# or
jupyter lab <path>
```
