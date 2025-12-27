---
name: job-watcher
description: Monitor and manage VASP jobs. Use when the user wants to check job status, diagnose failures, or set up job monitoring.
---

# VASP Job Watcher Skill

This skill helps monitor VASP calculations, diagnose issues, and suggest fixes.

## Checking Job Status

When asked to check a VASP job:

1. **Look for output files**:
   ```bash
   ls -la OUTCAR OSZICAR CONTCAR vasprun.xml 2>/dev/null
   ```

2. **Check if running**:
   ```bash
   # For SLURM
   squeue -u $USER
   # For PBS
   qstat -u $USER
   ```

3. **Check convergence**:
   ```python
   # Read OSZICAR for electronic convergence
   with open('OSZICAR') as f:
       lines = f.readlines()

   # Check for ionic convergence
   from ase.io import read
   try:
       atoms = read('OUTCAR')
       print("Calculation completed")
   except:
       print("Calculation still running or failed")
   ```

## Diagnosing Failures

When a job fails, check these in order:

### 1. Check OUTCAR for errors
```bash
grep -i "error\|warning\|fatal\|killed" OUTCAR | tail -20
```

### 2. Common failure patterns

| Pattern | Diagnosis | Fix |
|---------|-----------|-----|
| "ZBRENT: fatal error" | Electronic convergence | Increase NELM, adjust AMIX |
| "EDIFF not reached" | SCF not converged | Use ALGO=VeryFast, reduce EDIFF |
| "VERY BAD NEWS" | Force issues | Check structure, reduce POTIM |
| "internal error in subroutine" | Memory/array | Reduce NCORE, use LREAL=Auto |
| "POSCAR and POTCAR are incompatible" | Wrong POTCAR | Check element order |
| "killed" / "Segmentation fault" | Out of memory | Request more memory, reduce parallelization |

### 3. Check resources
```bash
# Memory usage (if job finished)
grep "maximum memory used" OUTCAR

# Time per ionic step
grep "LOOP+" OUTCAR | tail -5
```

## Auto-Fix Workflow

When fixing a failed job:

1. **Diagnose the issue**:
   ```python
   from vasp import Vasp
   calc = Vasp.read('.')  # Read from current directory

   # Check what went wrong
   with open('OUTCAR') as f:
       content = f.read()

   if 'ZBRENT' in content:
       fix = 'electronic_convergence'
   elif 'memory' in content.lower():
       fix = 'memory'
   # etc.
   ```

2. **Apply fix**:
   ```python
   # Electronic convergence fix
   if fix == 'electronic_convergence':
       calc.set(
           algo='VeryFast',
           nelm=200,
           amix=0.1,
           bmix=0.01,
       )

   # Memory fix
   elif fix == 'memory':
       calc.set(
           lreal='Auto',
           ncore=4,  # Reduce parallelization
       )
   ```

3. **Restart from CONTCAR**:
   ```python
   from ase.io import read

   # Read structure from CONTCAR if available
   if os.path.exists('CONTCAR') and os.path.getsize('CONTCAR') > 0:
       atoms = read('CONTCAR')
       calc.atoms = atoms
       calc.set(istart=1, icharg=1)  # Continue from wavefunctions
   ```

4. **Resubmit**:
   ```python
   # Restart calculation
   energy = calc.potential_energy
   ```

## Continuous Monitoring Pattern

For watching jobs over time:

```python
import time
from pathlib import Path

def watch_job(directory, check_interval=60):
    """Watch a VASP job and report status."""
    path = Path(directory)

    while True:
        # Check if job is done
        outcar = path / 'OUTCAR'
        if outcar.exists():
            content = outcar.read_text()

            if 'General timing and accounting' in content:
                print("✓ Job completed successfully")
                return 'completed'

            if any(err in content for err in ['error', 'ZBRENT', 'fatal']):
                print("✗ Job failed - diagnosing...")
                return 'failed'

        # Check if still in queue
        # (add scheduler-specific check here)

        print(f"Job running... (checked at {time.strftime('%H:%M:%S')})")
        time.sleep(check_interval)
```

## Job Recovery Script

Create a recovery script for common issues:

```python
#!/usr/bin/env python
"""recover_job.py - Attempt to recover a failed VASP job."""

import os
import sys
from pathlib import Path

def diagnose(outcar_path):
    """Diagnose the failure from OUTCAR."""
    content = Path(outcar_path).read_text()

    issues = []
    fixes = {}

    if 'ZBRENT' in content or 'EDIFF is reached' not in content:
        issues.append('Electronic convergence failed')
        fixes['nelm'] = 300
        fixes['algo'] = 'VeryFast'
        fixes['amix'] = 0.1

    if 'VERY BAD NEWS' in content:
        issues.append('Force calculation issues')
        fixes['potim'] = 0.1

    if 'memory' in content.lower() or 'killed' in content.lower():
        issues.append('Memory/resource issues')
        fixes['ncore'] = 4
        fixes['lreal'] = 'Auto'

    return issues, fixes

def recover(directory):
    """Attempt to recover the job."""
    outcar = Path(directory) / 'OUTCAR'
    contcar = Path(directory) / 'CONTCAR'

    issues, fixes = diagnose(outcar)

    print("Diagnosed issues:")
    for issue in issues:
        print(f"  - {issue}")

    print("\nSuggested fixes:")
    for key, value in fixes.items():
        print(f"  {key} = {value}")

    # Could auto-apply fixes here
    return fixes

if __name__ == '__main__':
    directory = sys.argv[1] if len(sys.argv) > 1 else '.'
    recover(directory)
```

## Asking Claude to Watch Jobs

When a user asks you to watch a job:

1. Ask for the job directory
2. Check current status
3. If failed, diagnose and suggest fixes
4. Offer to apply fixes and restart
5. Set up a monitoring loop if requested

Example interaction:
```
User: Watch my VASP job in /scratch/my_calc

Claude: I'll check on your VASP job...
[Reads OUTCAR, OSZICAR]

The job appears to have failed due to electronic convergence issues.
I found "ZBRENT: fatal error" in the OUTCAR.

Suggested fix:
- Increase NELM from 60 to 200
- Change ALGO from 'Normal' to 'VeryFast'
- Add AMIX=0.1, BMIX=0.01

Would you like me to apply these fixes and restart from CONTCAR?
```
