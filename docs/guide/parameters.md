# Parameters

The parameter system provides presets for common calculation types.

## Using Presets

```python
from vasp.parameters import get_preset, get_vdw_params, get_ldau_params

# Static calculation preset
static = get_preset('static')

# Relaxation preset
relax = get_preset('relax')

# VdW correction
vdw = get_vdw_params('d3bj')

# DFT+U
ldau = get_ldau_params(['Fe', 'O'], {'Fe': HubbardU(u=4.0)})
```

## Available Presets

### Calculation Types

| Preset | Description |
|--------|-------------|
| `static` | Single-point energy |
| `relax` | Ionic relaxation |
| `relax-cell` | Full cell relaxation |
| `md-nvt` | NVT molecular dynamics |
| `md-npt` | NPT molecular dynamics |

### Van der Waals

| Method | Description |
|--------|-------------|
| `d2` | Grimme D2 |
| `d3` | Grimme D3 (zero damping) |
| `d3bj` | Grimme D3-BJ (recommended) |
| `ts` | Tkatchenko-Scheffler |
| `ts-scs` | TS with self-consistent screening |

### Hybrid Functionals

| Method | Description |
|--------|-------------|
| `hse06` | HSE06 (recommended for solids) |
| `hse03` | HSE03 |
| `pbe0` | PBE0 |
| `b3lyp` | B3LYP |

### DFT+U

```python
from vasp.parameters import get_ldau_params, HubbardU

ldau = get_ldau_params(
    symbols=['Fe', 'O'],
    u_values={
        'Fe': HubbardU(u=4.0, j=0.0),
    }
)
```

## Merging Parameters

```python
# Combine multiple presets
params = {}
params.update(get_preset('static'))
params.update(get_vdw_params('d3bj'))

calc = Vasp(atoms=atoms, **params)
```
