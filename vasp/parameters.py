"""VASP parameter presets and validation.

Provides convenient parameter sets for common calculation types:
- Van der Waals corrections (DFT-D2, DFT-D3, TS, VdW-DF)
- DFT+U for strongly correlated systems
- Hybrid functionals (HSE06, PBE0, B3LYP)
- Spin-orbit coupling
- Machine learning force fields (MLFF)
- Molecular dynamics
- Phonon calculations
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

# =============================================================================
# Van der Waals Corrections
# =============================================================================

VDW_METHODS: dict[str, dict[str, Any]] = {
    # DFT-D2 (Grimme)
    "d2": {"ivdw": 10},
    "dft-d2": {"ivdw": 10},
    "grimme-d2": {"ivdw": 10},
    # DFT-D3 zero-damping
    "d3": {"ivdw": 11},
    "dft-d3": {"ivdw": 11},
    "d3-zero": {"ivdw": 11},
    # DFT-D3 with Becke-Johnson damping
    "d3bj": {"ivdw": 12},
    "d3-bj": {"ivdw": 12},
    "dft-d3-bj": {"ivdw": 12},
    # Tkatchenko-Scheffler
    "ts": {"ivdw": 2},
    "tkatchenko-scheffler": {"ivdw": 2},
    # TS with iterative Hirshfeld
    "ts-hirshfeld": {"ivdw": 20},
    "ts-ih": {"ivdw": 20},
    # Self-consistent screening (TS+SCS)
    "ts-scs": {"ivdw": 21},
    # Many-body dispersion
    "mbd": {"ivdw": 202},
    "mbd-rsscs": {"ivdw": 202},
    # dDsC dispersion correction
    "ddsc": {"ivdw": 4},
    # VdW-DF functionals (require LUSE_VDW)
    "vdw-df": {"luse_vdw": True, "aggac": 0.0, "gga": "RE"},
    "vdw-df2": {"luse_vdw": True, "aggac": 0.0, "gga": "ML", "zab_vdw": -1.8867},
    "optpbe-vdw": {"luse_vdw": True, "aggac": 0.0, "gga": "OR"},
    "optb88-vdw": {
        "luse_vdw": True,
        "aggac": 0.0,
        "gga": "BO",
        "param1": 0.1833333333,
        "param2": 0.2200000000,
    },
    "optb86b-vdw": {
        "luse_vdw": True,
        "aggac": 0.0,
        "gga": "MK",
        "param1": 0.1234,
        "param2": 1.0000,
    },
    "rev-vdw-df2": {
        "luse_vdw": True,
        "aggac": 0.0,
        "gga": "MK",
        "param1": 0.1234,
        "param2": 0.711357,
        "zab_vdw": -1.8867,
    },
}


def get_vdw_params(method: str) -> dict[str, Any]:
    """Get VASP parameters for van der Waals correction method.

    Args:
        method: Name of vdW method (e.g., 'd3', 'd3bj', 'ts', 'vdw-df2')

    Returns:
        Dict of VASP parameters for the vdW correction.

    Raises:
        ValueError: If method is not recognized.

    Example:
        >>> params = get_vdw_params('d3bj')
        >>> calc = Vasp(..., **params)
    """
    key = method.lower().replace("_", "-")
    if key not in VDW_METHODS:
        available = ", ".join(sorted(VDW_METHODS.keys()))
        raise ValueError(
            f"Unknown vdW method '{method}'.\n"
            f"Available methods: {available}\n\n"
            f"Examples:\n"
            f"  get_vdw_params('d3bj')  # DFT-D3 with BJ damping (recommended)\n"
            f"  get_vdw_params('d2')    # DFT-D2\n"
            f"  get_vdw_params('ts')    # Tkatchenko-Scheffler\n\n"
            f"See: https://www.vasp.at/wiki/index.php/IVDW"
        )
    return VDW_METHODS[key].copy()


# =============================================================================
# DFT+U Parameters
# =============================================================================


@dataclass
class HubbardU:
    """Hubbard U parameters for DFT+U calculations.

    Attributes:
        u: Coulomb U parameter in eV.
        j: Exchange J parameter in eV.
        l_angular: Angular momentum (0=s, 1=p, 2=d, 3=f).
    """

    u: float
    j: float = 0.0
    l_angular: int = 2  # Default to d-orbitals (l_angular=2)


# Common U values from literature (Dudarev method: U_eff = U - J)
COMMON_U_VALUES = {
    # 3d transition metals (d-orbitals, l_angular=2)
    "Ti": HubbardU(u=3.0, j=0.0, l_angular=2),
    "V": HubbardU(u=3.1, j=0.0, l_angular=2),
    "Cr": HubbardU(u=3.5, j=0.0, l_angular=2),
    "Mn": HubbardU(u=3.9, j=0.0, l_angular=2),
    "Fe": HubbardU(u=4.0, j=0.0, l_angular=2),
    "Co": HubbardU(u=3.3, j=0.0, l_angular=2),
    "Ni": HubbardU(u=6.4, j=0.0, l_angular=2),
    "Cu": HubbardU(u=4.0, j=0.0, l_angular=2),
    "Zn": HubbardU(u=7.5, j=0.0, l_angular=2),
    # 4d transition metals
    "Zr": HubbardU(u=2.0, j=0.0, l_angular=2),
    "Nb": HubbardU(u=2.0, j=0.0, l_angular=2),
    "Mo": HubbardU(u=2.0, j=0.0, l_angular=2),
    # Lanthanides (f-orbitals, l_angular=3)
    "Ce": HubbardU(u=5.0, j=0.0, l_angular=3),
    "Pr": HubbardU(u=5.0, j=0.0, l_angular=3),
    "Nd": HubbardU(u=5.0, j=0.0, l_angular=3),
    "Sm": HubbardU(u=5.0, j=0.0, l_angular=3),
    "Eu": HubbardU(u=5.0, j=0.0, l_angular=3),
    "Gd": HubbardU(u=5.0, j=0.0, l_angular=3),
    # Actinides (f-orbitals, l_angular=3)
    "U": HubbardU(u=4.0, j=0.0, l_angular=3),
    "Pu": HubbardU(u=4.0, j=0.0, l_angular=3),
    # p-block (for O 2p states in some oxides)
    "O": HubbardU(u=0.0, j=0.0, l_angular=1),
}


def get_ldau_params(
    symbols: list[str],
    u_values: dict[str, float | HubbardU] | None = None,
    ldautype: int = 2,
    ldauprint: int = 1,
    lmaxmix: int = 4,
) -> dict[str, Any]:
    """Generate DFT+U parameters for given atomic symbols.

    Uses the Dudarev (rotationally invariant) approach by default.

    Args:
        symbols: List of unique atomic symbols in POTCAR order.
        u_values: Dict mapping symbol to U value or HubbardU object.
            If None, uses COMMON_U_VALUES for known elements.
        ldautype: 1 = Liechtenstein, 2 = Dudarev (simplified).
        ldauprint: Verbosity (0, 1, or 2).
        lmaxmix: Max l for on-site density matrix mixing.

    Returns:
        Dict of VASP parameters for DFT+U.

    Example:
        >>> params = get_ldau_params(['Fe', 'O'], {'Fe': 4.0})
        >>> calc = Vasp(..., **params)
    """
    if u_values is None:
        u_values = {}

    ldauu = []
    ldauj = []
    ldaul = []

    for symbol in symbols:
        if symbol in u_values:
            val = u_values[symbol]
            if isinstance(val, HubbardU):
                ldauu.append(val.u)
                ldauj.append(val.j)
                ldaul.append(val.l_angular)
            else:
                # Simple float U value, use defaults
                hub = COMMON_U_VALUES.get(symbol, HubbardU(u=val, j=0.0, l_angular=2))
                ldauu.append(val)
                ldauj.append(hub.j)
                ldaul.append(hub.l_angular)
        elif symbol in COMMON_U_VALUES:
            hub = COMMON_U_VALUES[symbol]
            ldauu.append(hub.u)
            ldauj.append(hub.j)
            ldaul.append(hub.l_angular)
        else:
            # No U for this element
            ldauu.append(0.0)
            ldauj.append(0.0)
            ldaul.append(-1)  # -1 means no +U

    return {
        "ldau": True,
        "ldautype": ldautype,
        "ldaul": ldaul,
        "ldauu": ldauu,
        "ldauj": ldauj,
        "ldauprint": ldauprint,
        "lmaxmix": lmaxmix,
        "lasph": True,  # Non-spherical contributions important for +U
    }


# =============================================================================
# Hybrid Functionals
# =============================================================================

HYBRID_FUNCTIONALS = {
    # HSE06 (screened hybrid)
    "hse06": {
        "lhfcalc": True,
        "hfscreen": 0.2,
        "algo": "Damped",
        "time": 0.4,
        "precfock": "Fast",
    },
    "hse": {
        "lhfcalc": True,
        "hfscreen": 0.2,
        "algo": "Damped",
        "time": 0.4,
        "precfock": "Fast",
    },
    # HSE03 (original HSE)
    "hse03": {
        "lhfcalc": True,
        "hfscreen": 0.3,
        "algo": "Damped",
        "time": 0.4,
    },
    # PBE0
    "pbe0": {
        "lhfcalc": True,
        "aexx": 0.25,
        "algo": "Damped",
        "time": 0.4,
    },
    # B3LYP (not recommended for solids)
    "b3lyp": {
        "lhfcalc": True,
        "gga": "B3",
        "aexx": 0.2,
        "aggax": 0.72,
        "aggac": 0.81,
        "aldac": 0.19,
        "algo": "Damped",
        "time": 0.4,
    },
    # Range-separated hybrids
    "lc-wpbe": {
        "lhfcalc": True,
        "hfscreen": 0.4,  # Long-range only
        "aexx": 1.0,
        "algo": "Damped",
    },
}


def get_hybrid_params(
    functional: str,
    nkpts_factor: int = 1,
) -> dict[str, Any]:
    """Get VASP parameters for hybrid functional calculation.

    Args:
        functional: Name of hybrid functional (e.g., 'hse06', 'pbe0').
        nkpts_factor: Reduce k-points by this factor for HF part
            (useful for large cells).

    Returns:
        Dict of VASP parameters for hybrid calculation.

    Example:
        >>> params = get_hybrid_params('hse06')
        >>> calc = Vasp(..., xc='PBE', **params)
    """
    key = functional.lower()
    if key not in HYBRID_FUNCTIONALS:
        available = ", ".join(sorted(HYBRID_FUNCTIONALS.keys()))
        raise ValueError(
            f"Unknown hybrid functional '{functional}'.\n"
            f"Available functionals: {available}\n\n"
            f"Examples:\n"
            f"  get_hybrid_params('hse06')  # HSE06 (recommended for semiconductors)\n"
            f"  get_hybrid_params('pbe0')   # PBE0 (25% exact exchange)\n"
            f"  get_hybrid_params('b3lyp')  # B3LYP\n\n"
            f"Note: Hybrid calculations are computationally expensive!\n"
            f"See: https://www.vasp.at/wiki/index.php/Hybrid_functionals"
        )

    params = HYBRID_FUNCTIONALS[key].copy()

    if nkpts_factor > 1:
        params["nkred"] = nkpts_factor

    return params


# =============================================================================
# Spin-Orbit Coupling
# =============================================================================


def get_soc_params(
    saxis: tuple[float, float, float] = (0, 0, 1),
    lorbmom: bool = True,
    gga_compat: bool = True,
) -> dict[str, Any]:
    """Get parameters for spin-orbit coupling calculation.

    Note: Requires vasp_ncl binary compiled without -DNGZhalf.

    Args:
        saxis: Spin quantization axis (default: z-axis).
        lorbmom: Calculate orbital moments.
        gga_compat: Use GGA_COMPAT for gradient correction.

    Returns:
        Dict of VASP parameters for SOC.

    Example:
        >>> params = get_soc_params()
        >>> calc = Vasp(..., **params)  # Use vasp_ncl!
    """
    return {
        "lsorbit": True,
        "lnoncollinear": True,  # Automatically set by LSORBIT
        "saxis": list(saxis),
        "lorbmom": lorbmom,
        "gga_compat": gga_compat,
        "lmaxmix": 4,  # For d-electrons
    }


# =============================================================================
# Machine Learning Force Fields (VASP 6.3+)
# =============================================================================


@dataclass
class MLFFConfig:
    """Configuration for machine learning force field.

    Attributes:
        mode: MLFF mode (train, select, run, refit).
        rcut1: Cutoff for 2-body descriptors (Å).
        rcut2: Cutoff for 3-body descriptors (Å).
        mb: Max 2-body basis functions.
        mb3: Max 3-body basis functions.
        wtifor: Weight for forces in training.
        wtoten: Weight for energy in training.
        wtsif: Weight for stress in training.
        lmlff: Enable MLFF.
    """

    mode: str = "train"
    rcut1: float = 8.0
    rcut2: float = 4.0
    mb: int = 8000
    mb3: int = 8000
    wtifor: float = 10.0
    wtoten: float = 1.0
    wtsif: float = 1.0
    lmlff: bool = True


MLFF_MODES = {
    "train": "run",  # On-the-fly training
    "select": "select",  # Select training data
    "run": "run",  # Use existing FF
    "refit": "refit",  # Refit for fast prediction
}


def get_mlff_params(
    mode: str = "train", rcut1: float = 8.0, rcut2: float = 4.0, **kwargs
) -> dict[str, Any]:
    """Get parameters for machine learning force field.

    Requires VASP 6.3+ with MLFF support.

    Args:
        mode: 'train' (on-the-fly), 'run' (prediction), 'refit' (optimize).
        rcut1: Two-body cutoff radius in Å.
        rcut2: Three-body cutoff radius in Å.
        **kwargs: Additional ML_ parameters.

    Returns:
        Dict of VASP MLFF parameters.

    Example:
        >>> # On-the-fly training during MD
        >>> params = get_mlff_params('train')
        >>> calc = Vasp(..., ibrion=0, **params)
    """
    mode_key = mode.lower()
    if mode_key not in MLFF_MODES:
        raise ValueError(
            f"Unknown MLFF mode '{mode}'.\n"
            f"Valid modes: train, select, run, refit\n\n"
            f"Examples:\n"
            f"  get_mlff_params('train')  # Train new MLFF during MD\n"
            f"  get_mlff_params('run')    # Use existing MLFF for MD\n"
            f"  get_mlff_params('select') # Select structures for training\n\n"
            f"See: https://www.vasp.at/wiki/index.php/Category:Machine_learning"
        )

    params = {
        "ml_lmlff": True,
        "ml_mode": MLFF_MODES[mode_key],
        "ml_rcut1": rcut1,
        "ml_rcut2": rcut2,
    }

    # Add training weights for training mode
    if mode_key == "train":
        params.update(
            {
                "ml_wtifor": kwargs.get("wtifor", 10.0),
                "ml_wtoten": kwargs.get("wtoten", 1.0),
                "ml_wtsif": kwargs.get("wtsif", 1.0),
            }
        )

    # Add any additional ML_ parameters
    for key, val in kwargs.items():
        if key.startswith("ml_") or key.upper().startswith("ML_"):
            params[key.lower()] = val

    return params


# =============================================================================
# Molecular Dynamics
# =============================================================================

MD_ENSEMBLES = {
    "nve": {"mdalgo": 0, "smass": -3},
    "nvt-nose": {"mdalgo": 2, "smass": 0},  # Nose-Hoover
    "nvt-langevin": {"mdalgo": 3},
    "npt-parrinello": {"mdalgo": 3, "isif": 3},  # Parrinello-Rahman
}


def get_md_params(
    ensemble: str = "nvt-nose",
    temperature: float = 300.0,
    timestep: float = 1.0,  # fs
    nsteps: int = 1000,
    temperature_final: float | None = None,
) -> dict[str, Any]:
    """Get parameters for molecular dynamics.

    Args:
        ensemble: 'nve', 'nvt-nose', 'nvt-langevin', 'npt-parrinello'.
        temperature: Temperature in K.
        timestep: Time step in fs.
        nsteps: Number of MD steps.
        temperature_final: Final temperature for temperature ramp.

    Returns:
        Dict of VASP MD parameters.

    Example:
        >>> params = get_md_params('nvt-nose', temperature=500, nsteps=5000)
        >>> calc = Vasp(..., **params)
    """
    key = ensemble.lower()
    if key not in MD_ENSEMBLES:
        available = ", ".join(sorted(MD_ENSEMBLES.keys()))
        raise ValueError(
            f"Unknown MD ensemble '{ensemble}'.\n"
            f"Available ensembles: {available}\n\n"
            f"Examples:\n"
            f"  get_md_params('nve')          # Constant energy (microcanonical)\n"
            f"  get_md_params('nvt-nose')     # Nosé-Hoover thermostat\n"
            f"  get_md_params('nvt-langevin') # Langevin thermostat\n"
            f"  get_md_params('npt')          # Constant pressure\n\n"
            f"See: https://www.vasp.at/wiki/index.php/Molecular_dynamics"
        )

    params = {
        "ibrion": 0,  # MD
        "nsw": nsteps,
        "potim": timestep,
        "tebeg": temperature,
        "teend": temperature_final if temperature_final is not None else temperature,
        "nblock": 1,
        "kblock": nsteps,
        "isym": 0,  # No symmetry in MD
    }
    params.update(MD_ENSEMBLES[key])

    return params


# =============================================================================
# Phonon Calculations
# =============================================================================


def get_phonon_params(
    method: str = "dfpt",
    supercell: tuple[int, int, int] | None = None,
) -> dict[str, Any]:
    """Get parameters for phonon calculations.

    Args:
        method: 'dfpt' (IBRION=8), 'finite-diff' (IBRION=5/6), or 'phonopy'.
        supercell: Supercell size for finite differences (ignored for DFPT).

    Returns:
        Dict of VASP phonon parameters.

    Example:
        >>> params = get_phonon_params('dfpt')
        >>> calc = Vasp(..., **params)
    """
    method = method.lower()

    if method == "dfpt":
        return {
            "ibrion": 8,
            "lepsilon": True,  # Born charges and dielectric
            "nfree": 2,
            "addgrid": True,
        }
    elif method in ("finite-diff", "fd"):
        params = {
            "ibrion": 6,  # Central differences
            "nfree": 2,
            "potim": 0.015,  # Displacement in Å
            "isym": 0,
        }
        return params
    elif method == "phonopy":
        # Settings for Phonopy: just need forces
        return {
            "ibrion": -1,  # No VASP relaxation
            "nsw": 0,
            "isym": 0,
        }
    else:
        raise ValueError(
            f"Unknown phonon method '{method}'.\n"
            f"Valid methods: dfpt, finite-diff, phonopy\n\n"
            f"Examples:\n"
            f"  get_phonon_params('dfpt')       # Density functional perturbation theory\n"
            f"  get_phonon_params('finite-diff') # Finite differences (IBRION=5/6)\n"
            f"  get_phonon_params('phonopy')    # For use with Phonopy package\n\n"
            f"See: https://www.vasp.at/wiki/index.php/Phonons"
        )


# =============================================================================
# Optical Properties / GW / BSE
# =============================================================================


def get_optical_params(
    nbands: int | None = None,
    nedos: int = 2000,
) -> dict[str, Any]:
    """Get parameters for optical properties (frequency-dependent dielectric).

    Args:
        nbands: Number of bands (increase for accurate optics).
        nedos: DOS grid points.

    Returns:
        Dict of VASP parameters for optical calculation.
    """
    params = {
        "loptics": True,
        "cshift": 0.1,
        "nedos": nedos,
    }
    if nbands:
        params["nbands"] = nbands
    return params


def get_gw_params(
    algo: str = "gw0",
    nbands: int | None = None,
    nomega: int = 50,
    encutgw: float | None = None,
) -> dict[str, Any]:
    """Get parameters for GW quasiparticle calculation.

    Args:
        algo: 'gw0', 'gw', 'scgw0', 'scgw'.
        nbands: Number of bands for GW.
        nomega: Frequency grid points.
        encutgw: Response function cutoff.

    Returns:
        Dict of VASP GW parameters.

    Example:
        >>> params = get_gw_params('gw0', nbands=200)
        >>> calc = Vasp(..., **params)
    """
    algo_map = {
        "gw0": "GW0",
        "gw": "GW",
        "scgw0": "scGW0",
        "scgw": "scGW",
        "evgw0": "EVGW0",
        "evgw": "EVGW",
    }

    key = algo.lower()
    if key not in algo_map:
        available = ", ".join(sorted(algo_map.keys()))
        raise ValueError(
            f"Unknown GW algorithm '{algo}'.\n"
            f"Available algorithms: {available}\n\n"
            f"Examples:\n"
            f"  get_gw_params('gw0')   # Single-shot GW (G0W0)\n"
            f"  get_gw_params('scgw0') # Self-consistent GW0\n"
            f"  get_gw_params('evgw')  # Eigenvalue-only self-consistent GW\n\n"
            f"Note: GW calculations are very expensive!\n"
            f"See: https://www.vasp.at/wiki/index.php/GW_calculations"
        )

    params = {
        "algo": algo_map[key],
        "nomega": nomega,
        "loptics": True,
    }

    if nbands:
        params["nbands"] = nbands
    if encutgw:
        params["encutgw"] = encutgw

    return params


def get_bse_params(
    nbands: int | None = None,
    nbandso: int | None = None,
    nbandsv: int | None = None,
) -> dict[str, Any]:
    """Get parameters for BSE optical calculation.

    Args:
        nbands: Total number of bands.
        nbandso: Number of occupied bands in BSE.
        nbandsv: Number of unoccupied bands in BSE.

    Returns:
        Dict of VASP BSE parameters.
    """
    params = {
        "algo": "BSE",
        "antires": 0,  # Tamm-Dancoff approximation
    }

    if nbands:
        params["nbands"] = nbands
    if nbandso:
        params["nbandso"] = nbandso
    if nbandsv:
        params["nbandsv"] = nbandsv

    return params


# =============================================================================
# Calculation Presets
# =============================================================================

CALCULATION_PRESETS: dict[str, dict[str, Any]] = {
    "static": {
        "nsw": 0,
        "ibrion": -1,
    },
    "relax": {
        "nsw": 100,
        "ibrion": 2,
        "isif": 2,
        "ediffg": -0.02,
    },
    "relax-cell": {
        "nsw": 100,
        "ibrion": 2,
        "isif": 3,
        "ediffg": -0.02,
    },
    "band-structure": {
        "nsw": 0,
        "ibrion": -1,
        "icharg": 11,
        "lorbit": 11,
    },
    "dos": {
        "nsw": 0,
        "ibrion": -1,
        "icharg": 11,
        "lorbit": 11,
        "nedos": 2001,
        "ismear": -5,
    },
}


def get_preset(name: str) -> dict[str, Any]:
    """Get parameter preset for common calculation types.

    Args:
        name: Preset name ('static', 'relax', 'relax-cell', 'band-structure', 'dos').

    Returns:
        Dict of VASP parameters.
    """
    if name not in CALCULATION_PRESETS:
        available = ", ".join(sorted(CALCULATION_PRESETS.keys()))
        raise ValueError(
            f"Unknown preset '{name}'.\n"
            f"Available presets: {available}\n\n"
            f"Examples:\n"
            f"  get_preset('static')         # Single-point energy\n"
            f"  get_preset('relax')          # Geometry optimization (ions only)\n"
            f"  get_preset('relax-cell')     # Full relaxation (ions + cell)\n"
            f"  get_preset('dos')            # Density of states\n"
            f"  get_preset('band-structure') # Band structure calculation\n"
        )
    return CALCULATION_PRESETS[name].copy()
