# VASP Syntax Highlighting for VS Code

Syntax highlighting, tooltips, and auto-completion for VASP (Vienna Ab initio Simulation Package) input and output files.

## Features

- **Syntax Highlighting** for VASP files:
  - `INCAR` - Electronic and ionic parameters
  - `POSCAR` / `CONTCAR` - Atomic structure files
  - `KPOINTS` - K-point definitions
  - `OUTCAR` - Calculation output

- **Hover Tooltips** with parameter documentation:
  - Parameter descriptions
  - Type information
  - Default values
  - Allowed values
  - Links to official VASP Wiki

- **Auto-completion** for INCAR parameters with snippets

## Supported File Types

| File | Extensions/Names |
|------|-----------------|
| INCAR | `INCAR`, `INCAR.orig`, `*.incar` |
| POSCAR | `POSCAR`, `CONTCAR`, `*.poscar`, `*.vasp` |
| KPOINTS | `KPOINTS`, `*.kpoints` |
| OUTCAR | `OUTCAR`, `*.outcar` |

## Installation

### From Package

```bash
cd vscode-vasp-syntax
make install
```

### Development Mode

```bash
cd vscode-vasp-syntax
make install-dev
```

### Manual Installation

1. Build the extension:
   ```bash
   npm install
   npm run compile
   ```

2. Package as VSIX:
   ```bash
   npm install -g @vscode/vsce
   vsce package
   ```

3. Install in VS Code:
   ```bash
   code --install-extension vasp-syntax-1.0.0.vsix
   ```

## Usage

1. Open any VASP file (INCAR, POSCAR, etc.)
2. The file will automatically be recognized and highlighted
3. Hover over parameters to see documentation
4. Use Ctrl+Space for auto-completion in INCAR files

### Open VASP Wiki

Right-click on a parameter and select "Open VASP Wiki" or use the command palette:

1. Place cursor on a parameter name
2. Open Command Palette (Ctrl+Shift+P)
3. Run "VASP: Open Wiki Documentation"

## Supported Parameters

Over 100 VASP parameters are documented, organized by category:

- **Basic**: ENCUT, PREC, ALGO, EDIFF, NELM, etc.
- **K-points/Smearing**: ISMEAR, SIGMA, KSPACING, KGAMMA
- **Ionic**: NSW, IBRION, ISIF, EDIFFG, POTIM
- **Spin/Magnetism**: ISPIN, MAGMOM, NUPDOWN
- **DFT+U**: LDAU, LDAUTYPE, LDAUL, LDAUU, LDAUJ
- **Hybrid**: LHFCALC, HFSCREEN, AEXX
- **Van der Waals**: IVDW, LUSE_VDW, VDW_S6
- **Spin-Orbit**: LSORBIT, LNONCOLLINEAR, SAXIS
- **Output**: LWAVE, LCHARG, LORBIT, NEDOS
- **Parallelization**: NCORE, NPAR, KPAR
- **GW/BSE**: NOMEGA, ENCUTGW, NBANDSO
- **MLFF**: ML_LMLFF, ML_MODE, ML_RCUT1
- **NEB**: IMAGES, SPRING, LCLIMB

## Screenshots

### Syntax Highlighting

INCAR files are colorized by parameter category:

![INCAR highlighting](images/incar-highlight.png)

### Hover Tooltips

Hover over any parameter to see documentation:

![Hover tooltip](images/hover-tooltip.png)

### Auto-completion

Type to get suggestions with documentation:

![Completion](images/completion.png)

## Development

```bash
# Install dependencies
npm install

# Build
npm run compile

# Watch for changes
npm run watch

# Package
make package
```

## Contributing

Contributions are welcome! Please add new parameters to `src/extension.ts` in the `VASP_PARAMETERS` object.

## License

MIT

## Links

- [VASP Official Site](https://www.vasp.at/)
- [VASP Wiki](https://www.vasp.at/wiki/)
- [ASE-VASP Interface](https://github.com/jkitchin/vasp-ase)
