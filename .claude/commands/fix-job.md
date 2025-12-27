Diagnose and fix a failed VASP job, then restart it.

Arguments: $ARGUMENTS (job directory path, optional)

## Automatic Fix Workflow

1. **Identify the job directory**:
   - Use $ARGUMENTS if provided
   - Otherwise use current directory

2. **Read and diagnose OUTCAR**:
   ```bash
   grep -i "error\|fatal\|zbrent\|warning\|killed" OUTCAR | head -30
   ```

3. **Determine the failure type and apply fix**:

### Electronic Convergence (ZBRENT, EDIFF not reached)
```python
# Fixes to apply:
new_params = {
    'algo': 'VeryFast',  # or 'Damped' for hybrids
    'nelm': 200,
    'amix': 0.1,
    'bmix': 0.01,
    'ediff': 1e-5,  # Maybe loosen slightly
}
```

### Ionic Convergence (oscillating, not converging)
```python
new_params = {
    'ibrion': 1,     # Switch to quasi-Newton
    'potim': 0.1,    # Smaller step
    'ediffg': -0.03, # Loosen force criterion
}
```

### Memory Issues (killed, segfault)
```python
new_params = {
    'lreal': 'Auto',
    'ncore': 4,      # Reduce parallelization
    'kpar': 1,       # Reduce k-point parallelization
}
```

### Magnetic Convergence
```python
new_params = {
    'amix_mag': 0.8,
    'bmix_mag': 0.0001,
    'nelm': 200,
}
```

4. **Backup current files**:
   ```bash
   mkdir -p backup_$(date +%Y%m%d_%H%M%S)
   cp INCAR POSCAR OUTCAR OSZICAR backup_*/
   ```

5. **Restart from CONTCAR**:
   ```bash
   # If CONTCAR exists and has content
   if [ -s CONTCAR ]; then
       cp CONTCAR POSCAR
   fi
   ```

6. **Apply fixes to INCAR**:
   - Read current INCAR
   - Modify parameters
   - Write new INCAR

7. **Optional: Set continuation flags**:
   ```
   ISTART = 1   # Read WAVECAR
   ICHARG = 1   # Read CHGCAR
   ```

8. **Resubmit the job**:
   - For SLURM: `sbatch submit.sh`
   - For PBS: `qsub submit.sh`
   - For local: run VASP directly

## Example Session

```
User: Fix my failed job in /scratch/calc

Claude: Analyzing failed job in /scratch/calc...

Found error: ZBRENT: fatal error in BFGS - electronic convergence failed

Current settings:
  ALGO = Normal
  NELM = 60
  EDIFF = 1e-6

Applying fixes:
  ALGO = Normal → VeryFast
  NELM = 60 → 200
  Adding: AMIX = 0.1, BMIX = 0.01

CONTCAR found with 32 atoms - will restart from this geometry.

Backing up files to backup_20241227_143022/

[Creates new INCAR with fixes]
[Copies CONTCAR to POSCAR]

Ready to resubmit. Run: sbatch submit.sh
```
