---
name: troubleshoot
description: Troubleshoot VASP calculation issues. Use when the user encounters errors or unexpected results.
---

# VASP Troubleshooting Skill

## Common Error Categories

### 1. Electronic Convergence Issues

**Symptoms**: "ZBRENT: fatal error", "EDIFF not reached"

**Solutions**:
- Increase NELM (max electronic steps): `nelm=200`
- Use better mixing: `amix=0.1, bmix=0.01`
- For metals, use appropriate smearing: `ismear=1, sigma=0.1`
- Try different algorithm: `algo='Fast'` or `algo='VeryFast'`

### 2. Ionic Convergence Issues

**Symptoms**: Relaxation doesn't converge, oscillating forces

**Solutions**:
- Reduce POTIM: `potim=0.1` or lower
- Try different optimizer: `ibrion=1` (quasi-Newton)
- Check for symmetry issues
- Increase NSW if needed

### 3. Memory Issues

**Symptoms**: "REAL_OPTLAY: killed", segfaults

**Solutions**:
- Reduce NCORE/NPAR
- Use `lreal='Auto'` for large cells
- Reduce k-points if possible
- Request more memory in job script

### 4. Pseudopotential Issues

**Symptoms**: "POTCAR file not found"

**Solutions**:
- Check VASP_PP_PATH is set correctly
- Verify directory structure: `$VASP_PP_PATH/potpaw_PBE/{element}/POTCAR`
- Use correct element names (e.g., 'Fe' not 'Fe_sv' unless intended)

### 5. K-point Issues

**Symptoms**: Wrong energies, metallic behavior in insulators

**Solutions**:
- Increase k-point density
- For metals: use `ismear=1` or `ismear=2`
- For insulators/molecules: use `ismear=0` or `ismear=-5`
- Check k-point grid is compatible with symmetry

### 6. Magnetic Convergence

**Symptoms**: Oscillating magnetic moments, wrong ground state

**Solutions**:
- Set reasonable initial MAGMOM
- Use `amix_mag=0.8, bmix_mag=0.0001` for difficult cases
- Try different spin configurations
- For AFM, ensure correct supercell and initial moments

### 7. Hybrid Functional Issues

**Symptoms**: Very slow, doesn't converge

**Solutions**:
- Start from converged PBE calculation
- Use `algo='Damped'` or `algo='All'`
- Reduce k-points for hybrid part: `nkred=2`
- Use `precfock='Fast'` for screening

## Diagnostic Steps

1. **Check OUTCAR** for warnings:
   ```python
   with open('OUTCAR') as f:
       for line in f:
           if 'WARNING' in line or 'Error' in line:
               print(line)
   ```

2. **Check convergence**:
   ```python
   energies = []
   with open('OSZICAR') as f:
       for line in f:
           if 'E0' in line:
               energies.append(float(line.split()[4]))
   ```

3. **Check forces**:
   ```python
   forces = calc.forces
   max_force = np.abs(forces).max()
   print(f"Max force: {max_force:.4f} eV/Å")
   ```

## Quick Fixes

| Issue | Quick Fix |
|-------|-----------|
| Not converging | Increase NELM, adjust AMIX |
| Memory error | Reduce NCORE, use LREAL=Auto |
| Slow | Use ALGO=Fast, reduce k-points |
| Wrong magnetism | Check MAGMOM, use ISPIN=2 |
| Metallic insulator | Use correct ISMEAR |
