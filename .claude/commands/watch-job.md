Watch a VASP job and diagnose/fix issues if it fails.

Arguments: $ARGUMENTS (job directory path)

## Workflow

1. **Check job directory**:
   - If $ARGUMENTS provided, use that path
   - Otherwise, check current directory
   - Verify VASP files exist (INCAR, POSCAR, etc.)

2. **Assess current status**:
   ```bash
   # Check for output files
   ls -la OUTCAR OSZICAR CONTCAR vasprun.xml 2>/dev/null

   # Check last modification time
   stat OUTCAR 2>/dev/null | grep Modify
   ```

3. **Check if job completed**:
   ```bash
   grep "General timing and accounting" OUTCAR
   ```

4. **If job failed, diagnose**:
   ```bash
   # Look for errors
   grep -i "error\|warning\|fatal\|zbrent\|killed" OUTCAR | tail -20

   # Check convergence
   grep "EDIFF is reached" OUTCAR | wc -l
   tail -20 OSZICAR
   ```

5. **Suggest fixes based on diagnosis**:

   | Error | Fix |
   |-------|-----|
   | ZBRENT | ALGO=VeryFast, NELM=200, AMIX=0.1 |
   | Memory | LREAL=Auto, reduce NCORE |
   | Forces | Reduce POTIM, check structure |
   | K-points | Adjust ISMEAR for system type |

6. **Offer to apply fixes**:
   - Show the INCAR changes needed
   - Ask if user wants to restart from CONTCAR
   - Apply changes and restart if approved

7. **If job still running**:
   - Show current ionic/electronic step
   - Estimate time remaining
   - Report energy convergence trend

## Example Output

```
Checking VASP job in /scratch/my_calc...

Status: FAILED
Last modified: 2 hours ago

Diagnosis:
  ✗ Electronic convergence failed (ZBRENT error)
  ✓ Structure looks reasonable
  ✗ Only 45 of 60 NELM steps completed

Suggested fixes:
  1. Set ALGO = VeryFast (was: Normal)
  2. Set NELM = 200 (was: 60)
  3. Set AMIX = 0.1, BMIX = 0.01

CONTCAR exists with valid structure.

Would you like me to apply these fixes and restart?
```

## For Continuous Monitoring

If user wants ongoing monitoring:

```bash
# Check every 5 minutes
watch -n 300 'grep "LOOP+" OUTCAR | tail -1; grep "F=" OSZICAR | tail -1'
```

Or set up a simple monitor script.
