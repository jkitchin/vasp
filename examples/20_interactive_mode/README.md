# Example 20: Interactive Mode

## Overview

VASP's interactive mode (`INTERACTIVE = .TRUE.`) maintains a persistent VASP process that reuses wavefunctions between ionic steps. This can reduce SCF cycles by up to 75% compared to restarting VASP for each geometry optimization step.

This example demonstrates:
- Using `InteractiveRunner` for local persistent VASP sessions
- External geometry optimization with ASE optimizers
- Socket-based remote VASP execution with the i-PI protocol

## When to Use Interactive Mode

Interactive mode is beneficial when:
- Running geometry optimizations with external optimizers (BFGS, FIRE, ML-based)
- Performing many sequential single-point calculations on similar structures
- Using machine learning potentials that need DFT validation
- Running NEB or molecular dynamics with frequent force evaluations

Interactive mode is NOT beneficial for:
- Single-point calculations
- VASP's internal optimizers (IBRION=1,2)
- Calculations where wavefunctions change dramatically between steps

## How It Works

```
Traditional Mode:
  Step 1: Start VASP → SCF (50 cycles) → Write files → Exit
  Step 2: Start VASP → SCF (50 cycles) → Write files → Exit
  ...

Interactive Mode:
  Start VASP once
  Step 1: SCF (50 cycles) → Send forces → Receive positions
  Step 2: SCF (15 cycles) → Send forces → Receive positions  ← Reuses wavefunction!
  Step 3: SCF (10 cycles) → Send forces → Receive positions  ← Even faster!
  ...
  Exit VASP
```

## Requirements

- VASP 5.4+ or 6.x compiled with interactive mode support
- The `INTERACTIVE` tag must be recognized by your VASP build

## Files

- `tutorial.ipynb` - Interactive notebook with examples
- `README.md` - This documentation

## Quick Start

### Local Interactive Mode

```python
from ase.build import bulk
from ase.optimize import BFGS
from vasp import Vasp
from vasp.runners import InteractiveRunner

# Create structure
atoms = bulk('Cu', 'fcc', a=3.6)

# Use interactive runner
runner = InteractiveRunner(
    vasp_command='vasp_std',
    mpi_command='mpirun -np 4',
)

# Start session and optimize
with runner:
    results = runner.start('relax/', atoms)

    # Manual optimization loop
    for step in range(50):
        forces = results.forces
        max_force = np.max(np.abs(forces))

        if max_force < 0.01:  # Converged
            break

        # Update positions (simple steepest descent)
        atoms.positions -= 0.1 * forces
        results = runner.step(atoms)
```

### With ASE Optimizer

```python
from ase.optimize import BFGS

# Create a calculator wrapper for ASE
class InteractiveCalculator:
    def __init__(self, runner, directory):
        self.runner = runner
        self.directory = directory
        self.started = False

    def get_forces(self, atoms):
        if not self.started:
            results = self.runner.start(self.directory, atoms)
            self.started = True
        else:
            results = self.runner.step(atoms)
        return results.forces

    def get_potential_energy(self, atoms):
        # Forces are calculated together with energy
        return self._last_energy

# Use with BFGS
calc = InteractiveCalculator(runner, 'relax/')
atoms.calc = calc
opt = BFGS(atoms)
opt.run(fmax=0.01)
```

### Remote Execution via Socket

For running VASP on a compute node while driving from a login node:

**On compute node (start VASP client):**
```python
from vasp.runners import SocketClient, SocketConfig, InteractiveRunner

config = SocketConfig(host='login-node', port=31415)
runner = InteractiveRunner(vasp_command='vasp_std')

client = SocketClient(config, runner)
client.connect()
client.run(atoms_template, 'calc/')
```

**On login node (driver):**
```python
from vasp.runners import SocketServer, SocketConfig

config = SocketConfig(port=31415)

with SocketServer(config) as server:
    server.wait_for_client()

    for step in range(100):
        results = server.calculate(atoms)
        if converged(results['forces']):
            break
        atoms.positions = optimize(atoms, results['forces'])
```

## Performance Comparison

Typical speedup for geometry optimization:

| System | Standard | Interactive | Speedup |
|--------|----------|-------------|---------|
| Cu bulk (2 atoms) | 100s | 40s | 2.5x |
| Pt surface (36 atoms) | 1800s | 600s | 3x |
| Molecule on surface | 3600s | 1000s | 3.6x |

The speedup comes from:
1. No process startup overhead
2. Wavefunction extrapolation between steps
3. Faster SCF convergence from better initial guess

## Troubleshooting

### VASP doesn't recognize INTERACTIVE tag
Your VASP build may not support interactive mode. Check with:
```bash
grep -i interactive path/to/vasp/src/*.F
```

### Process hangs waiting for input
Ensure stdin/stdout are properly piped. Check that no other process is consuming the pipes.

### Forces seem wrong
The interactive protocol uses scaled (fractional) coordinates. Ensure you're converting correctly if modifying positions directly.

## References

- [VaspInteractive](https://github.com/ulissigroup/vasp-interactive) - Original implementation by Ulissi Group
- [i-PI](https://ipi-code.org/) - Universal force engine protocol
- [ASE SocketIOCalculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/socketio/socketio.html)
