import * as vscode from 'vscode';

// VASP parameter documentation with links to official VASP wiki
interface VaspParameter {
    description: string;
    type: string;
    default?: string;
    values?: string[];
    category: string;
    wikiPage: string;
}

const VASP_WIKI_BASE = 'https://www.vasp.at/wiki/index.php/';

const VASP_PARAMETERS: Record<string, VaspParameter> = {
    // ============================================================
    // Basic Parameters
    // ============================================================
    'SYSTEM': {
        description: 'Name of the system. Used as title in output files.',
        type: 'string',
        default: 'unknown system',
        category: 'Basic',
        wikiPage: 'SYSTEM'
    },
    'ENCUT': {
        description: 'Plane-wave energy cutoff in eV. Higher values give more accurate but slower calculations.',
        type: 'float',
        default: 'largest ENMAX in POTCAR',
        category: 'Basic',
        wikiPage: 'ENCUT'
    },
    'ENAUG': {
        description: 'Kinetic energy cutoff for the augmentation charges in eV.',
        type: 'float',
        default: 'EAUG from POTCAR',
        category: 'Basic',
        wikiPage: 'ENAUG'
    },
    'PREC': {
        description: 'Precision mode controlling FFT grids and other settings.',
        type: 'string',
        default: 'Normal',
        values: ['Low', 'Medium', 'High', 'Normal', 'Accurate', 'Single'],
        category: 'Basic',
        wikiPage: 'PREC'
    },
    'ALGO': {
        description: 'Electronic minimization algorithm.',
        type: 'string',
        default: 'Normal',
        values: ['Normal', 'VeryFast', 'Fast', 'Conjugate', 'All', 'Damped', 'Subrot', 'Eigenval', 'Exact', 'None', 'Nothing', 'CHI', 'GW0', 'GW', 'scGW0', 'scGW', 'EVGW0', 'EVGW', 'BSE'],
        category: 'Basic',
        wikiPage: 'ALGO'
    },
    'EDIFF': {
        description: 'Energy convergence criterion for electronic self-consistency in eV.',
        type: 'float',
        default: '1E-4',
        category: 'Basic',
        wikiPage: 'EDIFF'
    },
    'NELM': {
        description: 'Maximum number of electronic self-consistency steps.',
        type: 'integer',
        default: '60',
        category: 'Basic',
        wikiPage: 'NELM'
    },
    'NELMIN': {
        description: 'Minimum number of electronic self-consistency steps.',
        type: 'integer',
        default: '2',
        category: 'Basic',
        wikiPage: 'NELMIN'
    },
    'NELMDL': {
        description: 'Number of non-self-consistent steps at the beginning.',
        type: 'integer',
        default: '-5 (ISTART=0) or 0 (ISTART>0)',
        category: 'Basic',
        wikiPage: 'NELMDL'
    },
    'LREAL': {
        description: 'Determines whether projection operators are evaluated in real or reciprocal space.',
        type: 'string',
        default: '.FALSE.',
        values: ['.FALSE.', '.TRUE.', 'Auto', 'On', 'Off'],
        category: 'Basic',
        wikiPage: 'LREAL'
    },
    'ADDGRID': {
        description: 'Adds an additional support grid for augmentation charges.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Basic',
        wikiPage: 'ADDGRID'
    },
    'LASPH': {
        description: 'Include non-spherical contributions from gradient corrections inside PAW spheres.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Basic',
        wikiPage: 'LASPH'
    },
    'LPAW': {
        description: 'Use PAW (Projector Augmented Wave) method. Always .TRUE. for PAW potentials.',
        type: 'boolean',
        default: '.TRUE. (for PAW POTCAR)',
        category: 'Basic',
        wikiPage: 'LPAW'
    },
    'EAUG': {
        description: 'Augmentation charge cutoff from POTCAR (eV).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'EAUG'
    },
    'TITEL': {
        description: 'Title/name of the pseudopotential from POTCAR.',
        type: 'string',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'LULTRA': {
        description: 'Use ultrasoft pseudopotential.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'IUNSCR': {
        description: 'Unscreening type: 0=linear, 1=nonlinear, 2=none.',
        type: 'integer',
        values: ['0', '1', '2'],
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'RPACOR': {
        description: 'Partial core radius from POTCAR (a.u.).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'RCORE': {
        description: 'Outmost cutoff radius from POTCAR (a.u.).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'ENMAX': {
        description: 'Maximum recommended energy cutoff from POTCAR (eV).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'ENMIN': {
        description: 'Minimum recommended energy cutoff from POTCAR (eV).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'ICORE': {
        description: 'Core electron treatment from POTCAR.',
        type: 'integer',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'LCOR': {
        description: 'Correct augmentation charges.',
        type: 'boolean',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'VRHFIN': {
        description: 'Valence electron configuration from POTCAR.',
        type: 'string',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'GGA': {
        description: 'Type of GGA functional.',
        type: 'string',
        default: 'PE (PBE)',
        values: ['PE', 'RP', 'RE', 'PS', 'AM', 'MK', 'BO', 'OR', 'ML', 'B3', 'BF'],
        category: 'Basic',
        wikiPage: 'GGA'
    },
    'METAGGA': {
        description: 'Type of meta-GGA functional.',
        type: 'string',
        values: ['SCAN', 'TPSS', 'RTPSS', 'M06L', 'MBJ'],
        category: 'Basic',
        wikiPage: 'METAGGA'
    },
    'ISTART': {
        description: 'Determines how to start or restart a job.',
        type: 'integer',
        default: '0 (if no WAVECAR), 1 (if WAVECAR exists)',
        values: ['0', '1', '2'],
        category: 'Basic',
        wikiPage: 'ISTART'
    },
    'INIWAV': {
        description: 'Initial electronic wavefunction: 0=lowest eigenvalue of H, 1=random.',
        type: 'integer',
        default: '1',
        values: ['0', '1'],
        category: 'Basic',
        wikiPage: 'INIWAV'
    },
    'ICHARG': {
        description: 'Determines how the initial charge density is constructed.',
        type: 'integer',
        default: '2 (ISTART=0), 0 (ISTART>0)',
        values: ['0', '1', '2', '4', '10', '11', '12'],
        category: 'Basic',
        wikiPage: 'ICHARG'
    },

    // ============================================================
    // K-points and Smearing
    // ============================================================
    'ISMEAR': {
        description: 'Smearing method for partial occupancies. -5=tetrahedron, 0=Gaussian, 1-N=Methfessel-Paxton.',
        type: 'integer',
        default: '1',
        values: ['-5', '-4', '-3', '-2', '-1', '0', '1', '2'],
        category: 'Smearing',
        wikiPage: 'ISMEAR'
    },
    'SIGMA': {
        description: 'Width of smearing in eV.',
        type: 'float',
        default: '0.2',
        category: 'Smearing',
        wikiPage: 'SIGMA'
    },
    'KSPACING': {
        description: 'Distance between k-points in reciprocal space (1/Å). Automatically generates KPOINTS.',
        type: 'float',
        category: 'K-points',
        wikiPage: 'KSPACING'
    },
    'KGAMMA': {
        description: 'If .TRUE., the generated k-point grid is centered at the Gamma point.',
        type: 'boolean',
        default: '.TRUE.',
        category: 'K-points',
        wikiPage: 'KGAMMA'
    },

    // ============================================================
    // Ionic Relaxation
    // ============================================================
    'NSW': {
        description: 'Maximum number of ionic steps.',
        type: 'integer',
        default: '0',
        category: 'Ionic',
        wikiPage: 'NSW'
    },
    'IBRION': {
        description: 'Ionic relaxation algorithm. -1=no update, 0=MD, 1=quasi-Newton, 2=CG, 5-8=phonons.',
        type: 'integer',
        default: '-1 (NSW=0,1) or 0 (NSW>1)',
        values: ['-1', '0', '1', '2', '3', '5', '6', '7', '8', '44'],
        category: 'Ionic',
        wikiPage: 'IBRION'
    },
    'ISIF': {
        description: 'Stress tensor and what to relax. 0-1=atoms, 2=atoms+stress, 3=atoms+cell, 4-7=shape variants.',
        type: 'integer',
        default: '0 (IBRION=0) or 2 (IBRION>0)',
        values: ['0', '1', '2', '3', '4', '5', '6', '7'],
        category: 'Ionic',
        wikiPage: 'ISIF'
    },
    'EDIFFG': {
        description: 'Force convergence criterion. If negative, specifies max force in eV/Å.',
        type: 'float',
        default: 'EDIFF*10',
        category: 'Ionic',
        wikiPage: 'EDIFFG'
    },
    'POTIM': {
        description: 'Step width scaling for ionic movements (fs for MD, Å for relaxation).',
        type: 'float',
        default: '0.5',
        category: 'Ionic',
        wikiPage: 'POTIM'
    },
    'NFREE': {
        description: 'Number of displacements for finite difference phonon calculations.',
        type: 'integer',
        default: '2',
        values: ['1', '2', '4'],
        category: 'Ionic',
        wikiPage: 'NFREE'
    },
    'PSTRESS': {
        description: 'External pressure in kB. Used for enthalpy calculations (H = E + P*V).',
        type: 'float',
        default: '0',
        category: 'Ionic',
        wikiPage: 'PSTRESS'
    },

    // ============================================================
    // Spin and Magnetism
    // ============================================================
    'ISPIN': {
        description: 'Spin polarization: 1=non-spin-polarized, 2=spin-polarized.',
        type: 'integer',
        default: '1',
        values: ['1', '2'],
        category: 'Spin',
        wikiPage: 'ISPIN'
    },
    'MAGMOM': {
        description: 'Initial magnetic moment for each atom. For spin-polarized calculations.',
        type: 'float array',
        default: 'NIONS*1.0',
        category: 'Spin',
        wikiPage: 'MAGMOM'
    },
    'NUPDOWN': {
        description: 'Fix total magnetic moment to specified value (N_up - N_down).',
        type: 'integer',
        default: '-1 (not fixed)',
        category: 'Spin',
        wikiPage: 'NUPDOWN'
    },
    'LSORBIT': {
        description: 'Enable spin-orbit coupling. Requires non-collinear version of VASP.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Spin-Orbit',
        wikiPage: 'LSORBIT'
    },
    'LNONCOLLINEAR': {
        description: 'Enable non-collinear magnetism.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Spin-Orbit',
        wikiPage: 'LNONCOLLINEAR'
    },
    'SAXIS': {
        description: 'Quantization axis for spin (3D vector).',
        type: 'float array',
        default: '0 0 1',
        category: 'Spin-Orbit',
        wikiPage: 'SAXIS'
    },
    'LORBMOM': {
        description: 'Calculate and print orbital moments.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Spin-Orbit',
        wikiPage: 'LORBMOM'
    },

    // ============================================================
    // DFT+U
    // ============================================================
    'LDAU': {
        description: 'Enable DFT+U (LDA+U/GGA+U) calculations for strongly correlated systems.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'DFT+U',
        wikiPage: 'LDAU'
    },
    'LDAUTYPE': {
        description: 'Type of DFT+U approach: 1=Liechtenstein, 2=Dudarev (rotationally invariant).',
        type: 'integer',
        default: '2',
        values: ['1', '2', '4'],
        category: 'DFT+U',
        wikiPage: 'LDAUTYPE'
    },
    'LDAUL': {
        description: 'Angular momentum for each species (0=s, 1=p, 2=d, 3=f, -1=off).',
        type: 'integer array',
        category: 'DFT+U',
        wikiPage: 'LDAUL'
    },
    'LDAUU': {
        description: 'On-site Coulomb interaction U for each species (eV).',
        type: 'float array',
        category: 'DFT+U',
        wikiPage: 'LDAUU'
    },
    'LDAUJ': {
        description: 'On-site exchange interaction J for each species (eV).',
        type: 'float array',
        category: 'DFT+U',
        wikiPage: 'LDAUJ'
    },
    'LDAUPRINT': {
        description: 'Verbosity of DFT+U output: 0=silent, 1=occupancy matrix, 2=all.',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2'],
        category: 'DFT+U',
        wikiPage: 'LDAUPRINT'
    },
    'LMAXMIX': {
        description: 'Maximum l-quantum number for one-center PAW charge densities passed to mixer.',
        type: 'integer',
        default: '2 (s,p,d) or 4 (for f-elements)',
        values: ['2', '4', '6'],
        category: 'DFT+U',
        wikiPage: 'LMAXMIX'
    },

    // ============================================================
    // Hybrid Functionals
    // ============================================================
    'LHFCALC': {
        description: 'Enable hybrid functional calculation with Hartree-Fock exchange.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Hybrid',
        wikiPage: 'LHFCALC'
    },
    'HFSCREEN': {
        description: 'Screening length for HSE-type hybrid functionals (1/Å). 0.2 for HSE06.',
        type: 'float',
        default: '0',
        category: 'Hybrid',
        wikiPage: 'HFSCREEN'
    },
    'AEXX': {
        description: 'Fraction of exact exchange. 0.25 for PBE0/HSE06.',
        type: 'float',
        default: '0 (GGA) or 0.25 (hybrid)',
        category: 'Hybrid',
        wikiPage: 'AEXX'
    },
    'AGGAX': {
        description: 'Fraction of GGA exchange.',
        type: 'float',
        default: '1 - AEXX',
        category: 'Hybrid',
        wikiPage: 'AGGAX'
    },
    'AGGAC': {
        description: 'Fraction of GGA correlation.',
        type: 'float',
        default: '1.0',
        category: 'Hybrid',
        wikiPage: 'AGGAC'
    },
    'ALDAC': {
        description: 'Fraction of LDA correlation.',
        type: 'float',
        default: '1.0',
        category: 'Hybrid',
        wikiPage: 'ALDAC'
    },
    'TIME': {
        description: 'Time step for ALGO=Damped. Typically 0.4 for hybrids.',
        type: 'float',
        default: '0.4',
        category: 'Hybrid',
        wikiPage: 'TIME'
    },
    'PRECFOCK': {
        description: 'FFT grid for exact exchange: Fast, Normal, Accurate.',
        type: 'string',
        default: 'Normal',
        values: ['Low', 'Medium', 'Fast', 'Normal', 'Accurate'],
        category: 'Hybrid',
        wikiPage: 'PRECFOCK'
    },
    'NKRED': {
        description: 'Reduce k-points for HF by this factor in all directions.',
        type: 'integer',
        default: '1',
        category: 'Hybrid',
        wikiPage: 'NKRED'
    },

    // ============================================================
    // Van der Waals
    // ============================================================
    'IVDW': {
        description: 'Van der Waals correction method. 0=none, 10=DFT-D2, 11=DFT-D3, 12=DFT-D3(BJ), 2=TS, 202=MBD.',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2', '4', '10', '11', '12', '20', '21', '202', '263', '300'],
        category: 'Van der Waals',
        wikiPage: 'IVDW'
    },
    'LUSE_VDW': {
        description: 'Enable vdW-DF type non-local van der Waals functionals.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Van der Waals',
        wikiPage: 'LUSE_VDW'
    },
    'VDW_RADIUS': {
        description: 'Cutoff radius for vdW interactions (Å).',
        type: 'float',
        default: '50.0',
        category: 'Van der Waals',
        wikiPage: 'VDW_RADIUS'
    },
    'VDW_S6': {
        description: 'Global scaling parameter s6 for DFT-D2.',
        type: 'float',
        default: 'functional dependent',
        category: 'Van der Waals',
        wikiPage: 'VDW_S6'
    },
    'VDW_S8': {
        description: 'Scaling parameter s8 for DFT-D3.',
        type: 'float',
        default: 'functional dependent',
        category: 'Van der Waals',
        wikiPage: 'VDW_S8'
    },
    'VDW_A1': {
        description: 'Damping parameter a1 for DFT-D3(BJ).',
        type: 'float',
        category: 'Van der Waals',
        wikiPage: 'VDW_A1'
    },
    'VDW_A2': {
        description: 'Damping parameter a2 for DFT-D3(BJ).',
        type: 'float',
        category: 'Van der Waals',
        wikiPage: 'VDW_A2'
    },
    'PARAM1': {
        description: 'Parameter 1 for vdW-DF functionals.',
        type: 'float',
        category: 'Van der Waals',
        wikiPage: 'PARAM1'
    },
    'PARAM2': {
        description: 'Parameter 2 for vdW-DF functionals.',
        type: 'float',
        category: 'Van der Waals',
        wikiPage: 'PARAM2'
    },
    'ZAB_VDW': {
        description: 'vdW-DF kernel parameter. -1.8867 for vdW-DF2.',
        type: 'float',
        default: '-0.8491',
        category: 'Van der Waals',
        wikiPage: 'ZAB_VDW'
    },

    // ============================================================
    // Molecular Dynamics
    // ============================================================
    'MDALGO': {
        description: 'Molecular dynamics algorithm: 0=NVE, 1=Andersen, 2=Nose-Hoover, 3=Langevin.',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2', '3', '11', '13', '21'],
        category: 'Molecular Dynamics',
        wikiPage: 'MDALGO'
    },
    'SMASS': {
        description: 'Nose mass parameter. -3=NVE (microcanonical), >=0 for Nose-Hoover thermostat.',
        type: 'float',
        default: '-3',
        category: 'Molecular Dynamics',
        wikiPage: 'SMASS'
    },
    'TEBEG': {
        description: 'Initial temperature for MD (Kelvin).',
        type: 'float',
        default: '0',
        category: 'Molecular Dynamics',
        wikiPage: 'TEBEG'
    },
    'TEEND': {
        description: 'Final temperature for MD simulated annealing (Kelvin).',
        type: 'float',
        default: 'TEBEG',
        category: 'Molecular Dynamics',
        wikiPage: 'TEEND'
    },
    'NBLOCK': {
        description: 'After NBLOCK steps, pair correlation function and DOS are calculated.',
        type: 'integer',
        default: '1',
        category: 'Molecular Dynamics',
        wikiPage: 'NBLOCK'
    },
    'KBLOCK': {
        description: 'After KBLOCK*NBLOCK steps, averages are written to XDATCAR.',
        type: 'integer',
        default: 'NSW',
        category: 'Molecular Dynamics',
        wikiPage: 'KBLOCK'
    },
    'LANGEVIN_GAMMA': {
        description: 'Friction coefficient for Langevin thermostat (ps^-1) per species.',
        type: 'float array',
        category: 'Molecular Dynamics',
        wikiPage: 'LANGEVIN_GAMMA'
    },
    'PMASS': {
        description: 'Mass of the lattice degree of freedom for NPT dynamics.',
        type: 'float',
        category: 'Molecular Dynamics',
        wikiPage: 'PMASS'
    },

    // ============================================================
    // Output Control
    // ============================================================
    'LWAVE': {
        description: 'Write WAVECAR file (wavefunctions).',
        type: 'boolean',
        default: '.TRUE.',
        category: 'Output',
        wikiPage: 'LWAVE'
    },
    'LCHARG': {
        description: 'Write CHGCAR file (charge density).',
        type: 'boolean',
        default: '.TRUE.',
        category: 'Output',
        wikiPage: 'LCHARG'
    },
    'LVTOT': {
        description: 'Write LOCPOT file (local potential).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Output',
        wikiPage: 'LVTOT'
    },
    'LVHAR': {
        description: 'Write only Hartree potential to LOCPOT.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Output',
        wikiPage: 'LVHAR'
    },
    'LELF': {
        description: 'Write ELFCAR file (electron localization function).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Output',
        wikiPage: 'LELF'
    },
    'LAECHG': {
        description: 'Write all-electron charge density (AECCAR0, AECCAR2).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Output',
        wikiPage: 'LAECHG'
    },
    'LORBIT': {
        description: 'Write PROCAR/DOSCAR with projected DOS. 10-12 for atom-projected.',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2', '5', '10', '11', '12', '13', '14'],
        category: 'Output',
        wikiPage: 'LORBIT'
    },
    'NEDOS': {
        description: 'Number of grid points in DOS.',
        type: 'integer',
        default: '301',
        category: 'Output',
        wikiPage: 'NEDOS'
    },
    'EMIN': {
        description: 'Minimum energy for DOS/band evaluation (eV relative to E_Fermi).',
        type: 'float',
        default: 'lowest eigenvalue',
        category: 'Output',
        wikiPage: 'EMIN'
    },
    'EMAX': {
        description: 'Maximum energy for DOS/band evaluation (eV relative to E_Fermi).',
        type: 'float',
        default: 'highest eigenvalue',
        category: 'Output',
        wikiPage: 'EMAX'
    },
    'NWRITE': {
        description: 'Verbosity of OUTCAR. 0=minimal, 1=less, 2=normal, 3=verbose, 4=debug.',
        type: 'integer',
        default: '2',
        values: ['0', '1', '2', '3', '4'],
        category: 'Output',
        wikiPage: 'NWRITE'
    },
    'NBANDS': {
        description: 'Number of bands included in the calculation.',
        type: 'integer',
        default: 'depends on NELECT',
        category: 'Electronic',
        wikiPage: 'NBANDS'
    },
    'NELECT': {
        description: 'Total number of electrons in the calculation.',
        type: 'float',
        default: 'from POTCAR',
        category: 'Electronic',
        wikiPage: 'NELECT'
    },

    // ============================================================
    // Symmetry
    // ============================================================
    'ISYM': {
        description: 'Symmetry mode: -1=off, 0=on (no symmetry of k), 1=auto, 2=force, 3=no time-reversal.',
        type: 'integer',
        default: '1 or 2',
        values: ['-1', '0', '1', '2', '3'],
        category: 'Symmetry',
        wikiPage: 'ISYM'
    },
    'SYMPREC': {
        description: 'Precision for finding symmetry (Å).',
        type: 'float',
        default: '1E-5',
        category: 'Symmetry',
        wikiPage: 'SYMPREC'
    },

    // ============================================================
    // Parallelization
    // ============================================================
    'NCORE': {
        description: 'Number of cores per orbital. NCORE * NPAR = total cores.',
        type: 'integer',
        default: '1',
        category: 'Parallel',
        wikiPage: 'NCORE'
    },
    'NPAR': {
        description: 'Number of groups for band parallelization.',
        type: 'integer',
        default: 'total cores / NCORE',
        category: 'Parallel',
        wikiPage: 'NPAR'
    },
    'KPAR': {
        description: 'Number of k-point groups for k-point parallelization.',
        type: 'integer',
        default: '1',
        category: 'Parallel',
        wikiPage: 'KPAR'
    },
    'LPLANE': {
        description: 'Use plane-wave parallelization (usually .TRUE.).',
        type: 'boolean',
        default: '.TRUE.',
        category: 'Parallel',
        wikiPage: 'LPLANE'
    },
    'LSCALU': {
        description: 'Use parallel LU decomposition.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Parallel',
        wikiPage: 'LSCALU'
    },
    'NSIM': {
        description: 'Number of bands handled simultaneously in RMM-DIIS.',
        type: 'integer',
        default: '4',
        category: 'Parallel',
        wikiPage: 'NSIM'
    },

    // ============================================================
    // Dipole Corrections
    // ============================================================
    'LDIPOL': {
        description: 'Enable dipole corrections for surface calculations.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Dipole',
        wikiPage: 'LDIPOL'
    },
    'IDIPOL': {
        description: 'Direction of dipole correction: 1=x, 2=y, 3=z, 4=all directions.',
        type: 'integer',
        values: ['1', '2', '3', '4'],
        category: 'Dipole',
        wikiPage: 'IDIPOL'
    },
    'DIPOL': {
        description: 'Center of the cell for dipole correction (fractional coordinates).',
        type: 'float array',
        default: 'center of ionic charges',
        category: 'Dipole',
        wikiPage: 'DIPOL'
    },
    'EPSILON': {
        description: 'Dielectric constant of bulk for dipole correction.',
        type: 'float',
        default: '1.0',
        category: 'Dipole',
        wikiPage: 'EPSILON'
    },
    'EFIELD': {
        description: 'Applied electric field in eV/Å.',
        type: 'float',
        default: '0',
        category: 'Dipole',
        wikiPage: 'EFIELD'
    },

    // ============================================================
    // GW and BSE
    // ============================================================
    'NOMEGA': {
        description: 'Number of frequency points for GW/dielectric function.',
        type: 'integer',
        default: '0',
        category: 'GW/BSE',
        wikiPage: 'NOMEGA'
    },
    'ENCUTGW': {
        description: 'Energy cutoff for response function in GW (eV).',
        type: 'float',
        default: '2/3 * ENCUT',
        category: 'GW/BSE',
        wikiPage: 'ENCUTGW'
    },
    'NBANDSO': {
        description: 'Number of occupied bands in BSE calculation.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'NBANDSO'
    },
    'NBANDSV': {
        description: 'Number of virtual (unoccupied) bands in BSE calculation.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'NBANDSV'
    },
    'ANTIRES': {
        description: 'Antiresonant contributions in BSE: 0=with, 1=without (Tamm-Dancoff).',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2'],
        category: 'GW/BSE',
        wikiPage: 'ANTIRES'
    },
    'OMEGAMAX': {
        description: 'Maximum frequency for optical properties (eV).',
        type: 'float',
        default: '-30 to 30',
        category: 'GW/BSE',
        wikiPage: 'OMEGAMAX'
    },

    // ============================================================
    // Optical Properties
    // ============================================================
    'LOPTICS': {
        description: 'Calculate frequency-dependent dielectric function.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Optics',
        wikiPage: 'LOPTICS'
    },
    'CSHIFT': {
        description: 'Complex shift for dielectric function (eV). Broadening of spectra.',
        type: 'float',
        default: '0.1',
        category: 'Optics',
        wikiPage: 'CSHIFT'
    },
    'LEPSILON': {
        description: 'Calculate static dielectric tensor and Born effective charges.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Optics',
        wikiPage: 'LEPSILON'
    },
    'LRPA': {
        description: 'Include local field effects in RPA.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Optics',
        wikiPage: 'LRPA'
    },

    // ============================================================
    // Machine Learning Force Fields
    // ============================================================
    'ML_LMLFF': {
        description: 'Enable machine learning force field.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'MLFF',
        wikiPage: 'ML_LMLFF'
    },
    'ML_MODE': {
        description: 'MLFF operation mode: run, select, train, refit.',
        type: 'string',
        default: 'run',
        values: ['run', 'select', 'train', 'refit'],
        category: 'MLFF',
        wikiPage: 'ML_MODE'
    },
    'ML_ISTART': {
        description: 'When to start MLFF learning (ionic step).',
        type: 'integer',
        default: '0',
        category: 'MLFF',
        wikiPage: 'ML_ISTART'
    },
    'ML_RCUT1': {
        description: 'Cutoff radius for 2-body descriptors (Å).',
        type: 'float',
        default: '8.0',
        category: 'MLFF',
        wikiPage: 'ML_RCUT1'
    },
    'ML_RCUT2': {
        description: 'Cutoff radius for 3-body descriptors (Å).',
        type: 'float',
        default: '4.0',
        category: 'MLFF',
        wikiPage: 'ML_RCUT2'
    },
    'ML_MB': {
        description: 'Maximum number of 2-body basis functions.',
        type: 'integer',
        default: '8000',
        category: 'MLFF',
        wikiPage: 'ML_MB'
    },
    'ML_MB3': {
        description: 'Maximum number of 3-body basis functions.',
        type: 'integer',
        default: '8000',
        category: 'MLFF',
        wikiPage: 'ML_MB'
    },
    'ML_WTIFOR': {
        description: 'Weight for forces in MLFF training.',
        type: 'float',
        default: '10.0',
        category: 'MLFF',
        wikiPage: 'ML_WTIFOR'
    },
    'ML_WTOTEN': {
        description: 'Weight for total energy in MLFF training.',
        type: 'float',
        default: '1.0',
        category: 'MLFF',
        wikiPage: 'ML_WTOTEN'
    },
    'ML_WTSIF': {
        description: 'Weight for stress tensor in MLFF training.',
        type: 'float',
        default: '1.0',
        category: 'MLFF',
        wikiPage: 'ML_WTSIF'
    },

    // ============================================================
    // NEB (Nudged Elastic Band)
    // ============================================================
    'IMAGES': {
        description: 'Number of NEB images between initial and final states.',
        type: 'integer',
        category: 'NEB',
        wikiPage: 'IMAGES'
    },
    'SPRING': {
        description: 'Spring constant for NEB (eV/Å²).',
        type: 'float',
        default: '-5.0',
        category: 'NEB',
        wikiPage: 'SPRING'
    },
    'LCLIMB': {
        description: 'Enable climbing image NEB for better transition state.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NEB',
        wikiPage: 'LCLIMB'
    },
    'LTANGENTOLD': {
        description: 'Use old tangent definition for NEB.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NEB',
        wikiPage: 'LTANGENTOLD'
    },
    'LNEBCELL': {
        description: 'Allow cell shape to change in NEB.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NEB',
        wikiPage: 'LNEBCELL'
    },
    'ICHAIN': {
        description: 'Chain type for NEB: 0=normal, 1=reversal.',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2', '3'],
        category: 'NEB',
        wikiPage: 'ICHAIN'
    },

    // ============================================================
    // Phonons
    // ============================================================
    'DQ': {
        description: 'Displacement for finite difference phonon calculations (Å).',
        type: 'float',
        default: '0.015',
        category: 'Phonons',
        wikiPage: 'IBRION'
    },

    // ============================================================
    // Mixing Parameters
    // ============================================================
    'AMIX': {
        description: 'Linear mixing parameter for charge density.',
        type: 'float',
        default: '0.4',
        category: 'Mixing',
        wikiPage: 'AMIX'
    },
    'BMIX': {
        description: 'Cutoff wave vector for Kerker mixing scheme.',
        type: 'float',
        default: '1.0',
        category: 'Mixing',
        wikiPage: 'BMIX'
    },
    'AMIX_MAG': {
        description: 'Linear mixing parameter for magnetization density.',
        type: 'float',
        default: '1.6',
        category: 'Mixing',
        wikiPage: 'AMIX_MAG'
    },
    'AMIN': {
        description: 'Minimal mixing parameter.',
        type: 'float',
        default: '0.1',
        category: 'Mixing',
        wikiPage: 'AMIN'
    },
    'IMIX': {
        description: 'Mixing type: 0=no, 1=Kerker, 2=Tchebychev, 4=Broyden.',
        type: 'integer',
        default: '4',
        values: ['0', '1', '2', '4'],
        category: 'Mixing',
        wikiPage: 'IMIX'
    },
    'INIMIX': {
        description: 'Initial mixing type.',
        type: 'integer',
        default: '1',
        category: 'Mixing',
        wikiPage: 'INIMIX'
    },
    'MIXPRE': {
        description: 'Preconditioning for Broyden mixer.',
        type: 'integer',
        default: '1',
        category: 'Mixing',
        wikiPage: 'MIXPRE'
    },
    'MAXMIX': {
        description: 'Maximum number of steps stored in Broyden mixer.',
        type: 'integer',
        default: '-45',
        category: 'Mixing',
        wikiPage: 'MAXMIX'
    },
    'WC': {
        description: 'Weight factor for each step in Broyden mixing.',
        type: 'float',
        default: '1000',
        category: 'Mixing',
        wikiPage: 'WC'
    },

    // ============================================================
    // Algorithm Parameters
    // ============================================================
    'IALGO': {
        description: 'Algorithm for electronic minimization (deprecated, use ALGO).',
        type: 'integer',
        default: '38',
        category: 'Algorithm',
        wikiPage: 'IALGO'
    },
    'LDIAG': {
        description: 'Sub-space diagonalization in eigenvalue solver.',
        type: 'boolean',
        default: '.TRUE.',
        category: 'Algorithm',
        wikiPage: 'LDIAG'
    },
    'LSUBROT': {
        description: 'Optimize orbitals by subspace rotation.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Algorithm',
        wikiPage: 'LSUBROT'
    },
    'WEIMIN': {
        description: 'Maximum weight for eigenvalue minimization.',
        type: 'float',
        default: '0.001',
        category: 'Algorithm',
        wikiPage: 'WEIMIN'
    },
    'EBREAK': {
        description: 'Break condition for eigenvalue optimization.',
        type: 'float',
        category: 'Algorithm',
        wikiPage: 'EBREAK'
    },
    'DEPER': {
        description: 'Relative energy change for break condition.',
        type: 'float',
        default: '0.3',
        category: 'Algorithm',
        wikiPage: 'DEPER'
    },
    'ENINI': {
        description: 'Initial cutoff energy for wavefunctions.',
        type: 'float',
        category: 'Algorithm',
        wikiPage: 'ENINI'
    },

    // ============================================================
    // Grid Parameters
    // ============================================================
    'NGXF': {
        description: 'FFT grid points in x for augmentation charges.',
        type: 'integer',
        category: 'Grid',
        wikiPage: 'NGXF'
    },
    'NGYF': {
        description: 'FFT grid points in y for augmentation charges.',
        type: 'integer',
        category: 'Grid',
        wikiPage: 'NGYF'
    },
    'NGZF': {
        description: 'FFT grid points in z for augmentation charges.',
        type: 'integer',
        category: 'Grid',
        wikiPage: 'NGZF'
    },
    'NGX': {
        description: 'FFT grid points in x direction.',
        type: 'integer',
        category: 'Grid',
        wikiPage: 'NGX'
    },
    'NGY': {
        description: 'FFT grid points in y direction.',
        type: 'integer',
        category: 'Grid',
        wikiPage: 'NGY'
    },
    'NGZ': {
        description: 'FFT grid points in z direction.',
        type: 'integer',
        category: 'Grid',
        wikiPage: 'NGZ'
    },

    // ============================================================
    // Additional Common Parameters
    // ============================================================
    'APACO': {
        description: 'Distance for pair correlation function (Å).',
        type: 'float',
        default: '16.0',
        category: 'MD',
        wikiPage: 'APACO'
    },
    'NPACO': {
        description: 'Number of slots for pair correlation function.',
        type: 'integer',
        default: '256',
        category: 'MD',
        wikiPage: 'NPACO'
    },
    'NLSPLINE': {
        description: 'Spline interpolation for projectors.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Basic',
        wikiPage: 'NLSPLINE'
    },
    'LCOMPAT': {
        description: 'Compatibility mode with VASP 4.4.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Basic',
        wikiPage: 'LCOMPAT'
    },
    'GGA_COMPAT': {
        description: 'GGA compatibility with VASP 4.4-4.6.',
        type: 'boolean',
        default: '.TRUE.',
        category: 'Basic',
        wikiPage: 'GGA_COMPAT'
    },
    'LMAXPAW': {
        description: 'Maximum l for on-site density.',
        type: 'integer',
        default: '-100 (auto)',
        category: 'PAW',
        wikiPage: 'LMAXPAW'
    },
    'ROPT': {
        description: 'Real-space optimization parameters per species.',
        type: 'float array',
        category: 'Basic',
        wikiPage: 'ROPT'
    },
    'VOSKOWN': {
        description: 'Vosko-Wilk-Nusair interpolation: 0=Perdew, 1=VWN.',
        type: 'integer',
        default: '0',
        values: ['0', '1'],
        category: 'XC',
        wikiPage: 'VOSKOWN'
    },

    // ============================================================
    // POTCAR / Pseudopotential Parameters
    // ============================================================
    'POMASS': {
        description: 'Ionic mass for each species from POTCAR (in atomic mass units).',
        type: 'float array',
        category: 'POTCAR',
        wikiPage: 'POMASS'
    },
    'ZVAL': {
        description: 'Number of valence electrons for each species from POTCAR.',
        type: 'float array',
        category: 'POTCAR',
        wikiPage: 'ZVAL'
    },
    'RWIGS': {
        description: 'Wigner-Seitz radius for each species (Å). Used for LORBIT projected DOS.',
        type: 'float array',
        category: 'POTCAR',
        wikiPage: 'RWIGS'
    },
    'EATOM': {
        description: 'Atomic reference energy from POTCAR (eV).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'RMAX': {
        description: 'Maximum radius for radial grids in POTCAR (a.u.).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'RAUG': {
        description: 'Augmentation sphere radius from POTCAR (a.u.).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'RDEP': {
        description: 'Core radius for depletion charge from POTCAR (a.u.).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'RDEPT': {
        description: 'Core radius for depletion charge kinetic energy (a.u.).',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'QCUT': {
        description: 'Plane wave cutoff for POTCAR projectors.',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'QGAM': {
        description: 'Gamma for POTCAR optimization.',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'DEXC': {
        description: 'Exchange-correlation energy difference from POTCAR.',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'EREF': {
        description: 'Reference energy from POTCAR.',
        type: 'float',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },
    'LEXCH': {
        description: 'Exchange-correlation type in POTCAR (e.g., PE=PBE, CA=LDA).',
        type: 'string',
        category: 'POTCAR',
        wikiPage: 'POTCAR'
    },

    // ============================================================
    // Hybrid Functional Additional Parameters
    // ============================================================
    'HFSCREENC': {
        description: 'Screening for correlation in hybrid functionals.',
        type: 'float',
        default: '0',
        category: 'Hybrid',
        wikiPage: 'HFSCREEN'
    },
    'HFRCUT': {
        description: 'Real-space cutoff for HF exchange (Å).',
        type: 'float',
        default: '0',
        category: 'Hybrid',
        wikiPage: 'HFRCUT'
    },
    'HFALPHA': {
        description: 'Yukawa potential screening for exchange.',
        type: 'float',
        default: '0',
        category: 'Hybrid',
        wikiPage: 'HFALPHA'
    },
    'HFKIDENT': {
        description: 'Use k-point identity for HF.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Hybrid',
        wikiPage: 'HFKIDENT'
    },
    'LRHFCALC': {
        description: 'Long-range Hartree-Fock calculation.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Hybrid',
        wikiPage: 'LRHFCALC'
    },
    'LHARTREE': {
        description: 'Include Hartree energy in HF.',
        type: 'boolean',
        default: '.TRUE.',
        category: 'Hybrid',
        wikiPage: 'LHFCALC'
    },
    'ENCUTFOCK': {
        description: 'FFT grid cutoff for exact exchange (eV).',
        type: 'float',
        category: 'Hybrid',
        wikiPage: 'ENCUTFOCK'
    },
    'LMAXFOCK': {
        description: 'Maximum L for charge augmentation in HF.',
        type: 'integer',
        category: 'Hybrid',
        wikiPage: 'LMAXFOCK'
    },
    'LMAXFOCKAE': {
        description: 'Maximum L for AE charge augmentation in HF.',
        type: 'integer',
        category: 'Hybrid',
        wikiPage: 'LMAXFOCKAE'
    },
    'NMAXFOCKAE': {
        description: 'Maximum number of AE augmentation channels.',
        type: 'integer',
        category: 'Hybrid',
        wikiPage: 'NMAXFOCKAE'
    },
    'ALDAX': {
        description: 'Fraction of LDA exchange (usually 0 for hybrids).',
        type: 'float',
        default: '0',
        category: 'Hybrid',
        wikiPage: 'ALDAX'
    },
    'EXXOEP': {
        description: 'Exact exchange OEP method: 0=standard, 1=local, 2=mixed.',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2'],
        category: 'Hybrid',
        wikiPage: 'EXXOEP'
    },
    'LMODELHF': {
        description: 'Use model HF (Thomas-Fermi screening).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Hybrid',
        wikiPage: 'LMODELHF'
    },

    // ============================================================
    // GW and BSE Additional Parameters
    // ============================================================
    'NBANDSGW': {
        description: 'Number of bands for GW calculation.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'NBANDSGW'
    },
    'NBANDSGWLOW': {
        description: 'Lowest band for GW calculation.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'NBANDSGWLOW'
    },
    'ENCUTGWSOFT': {
        description: 'Soft cutoff for response function (eV).',
        type: 'float',
        category: 'GW/BSE',
        wikiPage: 'ENCUTGWSOFT'
    },
    'ENCUTLF': {
        description: 'Energy cutoff for local field effects (eV).',
        type: 'float',
        category: 'GW/BSE',
        wikiPage: 'ENCUTLF'
    },
    'ENCUT4O': {
        description: 'Cutoff for 4-orbital integrals (eV).',
        type: 'float',
        category: 'GW/BSE',
        wikiPage: 'ENCUT4O'
    },
    'NOMEGAR': {
        description: 'Number of frequency points for real-axis.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'NOMEGAR'
    },
    'OMEGAMIN': {
        description: 'Minimum frequency for dielectric function (eV).',
        type: 'float',
        category: 'GW/BSE',
        wikiPage: 'OMEGAMIN'
    },
    'OMEGAGRID': {
        description: 'Type of frequency grid for GW.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'OMEGAGRID'
    },
    'OMEGATL': {
        description: 'Maximum frequency for time-dependent calculations (eV).',
        type: 'float',
        category: 'GW/BSE',
        wikiPage: 'OMEGATL'
    },
    'OMEGAPAR': {
        description: 'Parallelization over frequency points.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'OMEGAPAR'
    },
    'LFERMIGW': {
        description: 'Update Fermi energy in GW.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LFERMIGW'
    },
    'LSPECTRAL': {
        description: 'Use spectral method for dielectric function.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LSPECTRAL'
    },
    'LSPECTRALGW': {
        description: 'Calculate spectral function in GW.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LSPECTRALGW'
    },
    'SELFENERGY': {
        description: 'Self-energy evaluation method.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'SELFENERGY'
    },
    'NATURALO': {
        description: 'Natural orbital occupations.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'NATURALO'
    },
    'LSINGLES': {
        description: 'Include single excitations.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LSINGLES'
    },
    'LTRIPLET': {
        description: 'Calculate triplet excitations in BSE.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LTRIPLET'
    },
    'LADDER': {
        description: 'Include ladder diagrams.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LADDER'
    },
    'LFXC': {
        description: 'Include exchange-correlation kernel.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LFXC'
    },
    'LFXCEPS': {
        description: 'Include fxc in dielectric function.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LFXCEPS'
    },
    'LFXHEG': {
        description: 'Use homogeneous electron gas fxc.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LFXHEG'
    },
    'IBSE': {
        description: 'BSE type: 0=direct, 1=iterative, 2=Haydock.',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2'],
        category: 'GW/BSE',
        wikiPage: 'IBSE'
    },
    'LRSRPA': {
        description: 'Long-range RPA correlation.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LRSRPA'
    },
    'LRSCOR': {
        description: 'Long-range correlation.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'LRSCOR'
    },
    'NELMHF': {
        description: 'Number of electronic steps for HF.',
        type: 'integer',
        default: '1',
        category: 'GW/BSE',
        wikiPage: 'NELMHF'
    },
    'SCISSOR': {
        description: 'Scissors operator shift (eV).',
        type: 'float',
        default: '0',
        category: 'GW/BSE',
        wikiPage: 'SCISSOR'
    },
    'L2ORDER': {
        description: 'Second order contribution in RPA.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'L2ORDER'
    },

    // ============================================================
    // Van der Waals Additional Parameters
    // ============================================================
    'BPARAM': {
        description: 'B parameter for DFT-D corrections.',
        type: 'float',
        category: 'Van der Waals',
        wikiPage: 'DFT-D2'
    },
    'CPARAM': {
        description: 'C parameter for DFT-D corrections.',
        type: 'float',
        category: 'Van der Waals',
        wikiPage: 'DFT-D2'
    },
    'LVDWSCS': {
        description: 'Self-consistent screening for vdW.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Van der Waals',
        wikiPage: 'LVDWSCS'
    },

    // ============================================================
    // Berry Phase and Electric Field
    // ============================================================
    'LBERRY': {
        description: 'Calculate Berry phase for polarization.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Berry Phase',
        wikiPage: 'LBERRY'
    },
    'IGPAR': {
        description: 'Direction for Berry phase: 1=a, 2=b, 3=c.',
        type: 'integer',
        values: ['1', '2', '3'],
        category: 'Berry Phase',
        wikiPage: 'IGPAR'
    },
    'NPPSTR': {
        description: 'Number of k-points along Berry phase direction.',
        type: 'integer',
        category: 'Berry Phase',
        wikiPage: 'NPPSTR'
    },
    'LCALCEPS': {
        description: 'Calculate dielectric tensor.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Berry Phase',
        wikiPage: 'LCALCEPS'
    },
    'LCALCPOL': {
        description: 'Calculate polarization.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Berry Phase',
        wikiPage: 'LCALCPOL'
    },
    'LPEAD': {
        description: 'Derivative of orbitals w.r.t. k for optical properties.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Berry Phase',
        wikiPage: 'LPEAD'
    },

    // ============================================================
    // Solvation and Implicit Solvent
    // ============================================================
    'LSOL': {
        description: 'Enable implicit solvation model (VASPsol).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Solvation',
        wikiPage: 'LSOL'
    },
    'EB_K': {
        description: 'Bulk dielectric constant of solvent.',
        type: 'float',
        default: '78.4',
        category: 'Solvation',
        wikiPage: 'EB_K'
    },
    'TAU': {
        description: 'Surface tension parameter for solvation (eV/Å²).',
        type: 'float',
        category: 'Solvation',
        wikiPage: 'TAU'
    },
    'LAMBDA': {
        description: 'Debye screening length for solvation.',
        type: 'float',
        category: 'Solvation',
        wikiPage: 'LAMBDA'
    },

    // ============================================================
    // Core Level and Spectroscopy
    // ============================================================
    'ICORELEVEL': {
        description: 'Core level shift calculation: 0=none, 1=initial, 2=final.',
        type: 'integer',
        default: '0',
        values: ['0', '1', '2'],
        category: 'Core Level',
        wikiPage: 'ICORELEVEL'
    },
    'CLNT': {
        description: 'Principal quantum number for core level.',
        type: 'integer',
        category: 'Core Level',
        wikiPage: 'CLNT'
    },
    'CLN': {
        description: 'Principal quantum number n for core hole.',
        type: 'integer',
        category: 'Core Level',
        wikiPage: 'CLN'
    },
    'CLL': {
        description: 'Angular momentum l for core hole.',
        type: 'integer',
        category: 'Core Level',
        wikiPage: 'CLL'
    },
    'CLZ': {
        description: 'Effective nuclear charge for core level.',
        type: 'float',
        category: 'Core Level',
        wikiPage: 'CLZ'
    },

    // ============================================================
    // NMR and Magnetic Response
    // ============================================================
    'LCHIMAG': {
        description: 'Calculate magnetic susceptibility.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NMR',
        wikiPage: 'LCHIMAG'
    },
    'LNMR_SYM_RED': {
        description: 'Reduce symmetry for NMR calculations.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NMR',
        wikiPage: 'LNMR_SYM_RED'
    },
    'ORBITALMAG': {
        description: 'Calculate orbital magnetization.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NMR',
        wikiPage: 'ORBITALMAG'
    },
    'LNABLA': {
        description: 'Use nabla operator for momentum matrix.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NMR',
        wikiPage: 'LNABLA'
    },
    'LEFG': {
        description: 'Calculate electric field gradient.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NMR',
        wikiPage: 'LEFG'
    },
    'LLRAUG': {
        description: 'Include augmentation for linear response.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NMR',
        wikiPage: 'LLRAUG'
    },
    'LRHOB': {
        description: 'Include B-field in density.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NMR',
        wikiPage: 'LRHOB'
    },

    // ============================================================
    // Spiral and Non-Collinear Magnetism
    // ============================================================
    'LSPIRAL': {
        description: 'Generalized Bloch spin spiral calculation.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Spin Spiral',
        wikiPage: 'LSPIRAL'
    },
    'QSPIRAL': {
        description: 'Spin spiral propagation vector.',
        type: 'float array',
        category: 'Spin Spiral',
        wikiPage: 'QSPIRAL'
    },
    'LZEROZ': {
        description: 'Force zero z-component of magnetization.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Spin Spiral',
        wikiPage: 'LZEROZ'
    },
    'LMAGBLOCH': {
        description: 'Calculate Bloch representation of magnetization.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Spin Spiral',
        wikiPage: 'LMAGBLOCH'
    },

    // ============================================================
    // Parallelization Additional Parameters
    // ============================================================
    'NKREDX': {
        description: 'Reduce k-points for HF in x direction.',
        type: 'integer',
        default: '1',
        category: 'Parallel',
        wikiPage: 'NKREDX'
    },
    'NKREDY': {
        description: 'Reduce k-points for HF in y direction.',
        type: 'integer',
        default: '1',
        category: 'Parallel',
        wikiPage: 'NKREDY'
    },
    'NKREDZ': {
        description: 'Reduce k-points for HF in z direction.',
        type: 'integer',
        default: '1',
        category: 'Parallel',
        wikiPage: 'NKREDZ'
    },
    'NKREDLFX': {
        description: 'Reduce k-points for local field effects.',
        type: 'integer',
        default: '1',
        category: 'Parallel',
        wikiPage: 'NKREDLFX'
    },

    // ============================================================
    // Advanced Algorithm Parameters
    // ============================================================
    'LBFGS': {
        description: 'Use L-BFGS optimizer (requires IBRION=3).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Algorithm',
        wikiPage: 'IBRION'
    },
    'LHF': {
        description: 'Enable Hartree-Fock (same as LHFCALC).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Hybrid',
        wikiPage: 'LHFCALC'
    },
    'LIBXC': {
        description: 'Use LibXC library for XC functionals.',
        type: 'integer',
        category: 'XC',
        wikiPage: 'LIBXC'
    },
    'ICHIBARE': {
        description: 'Bare susceptibility calculation type.',
        type: 'integer',
        category: 'Response',
        wikiPage: 'ICHIBARE'
    },
    'SMIX': {
        description: 'Mixing parameter for Kerker scheme.',
        type: 'float',
        default: '0.4',
        category: 'Mixing',
        wikiPage: 'SMIX'
    },
    'IWAVPR': {
        description: 'Wavefunction prediction: 0=none, 1=charge, 2=wave, 3=both.',
        type: 'integer',
        default: '2',
        values: ['0', '1', '2', '3'],
        category: 'Algorithm',
        wikiPage: 'IWAVPR'
    },
    'IRESTART': {
        description: 'Restart mode for MD.',
        type: 'integer',
        default: '0',
        category: 'MD',
        wikiPage: 'IRESTART'
    },
    'NREBOOT': {
        description: 'Number of reboots for MD.',
        type: 'integer',
        default: '0',
        category: 'MD',
        wikiPage: 'NREBOOT'
    },
    'SCALEE': {
        description: 'Scaling for kinetic energy in MD.',
        type: 'float',
        default: '1.0',
        category: 'MD',
        wikiPage: 'SCALEE'
    },
    'LMONO': {
        description: 'Monopole corrections.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Dipole',
        wikiPage: 'LMONO'
    },
    'LASYNC': {
        description: 'Asynchronous communication.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Parallel',
        wikiPage: 'LASYNC'
    },
    'LINTERFAST': {
        description: 'Fast interpolation.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Algorithm',
        wikiPage: 'LINTERFAST'
    },

    // ============================================================
    // DFPT and Linear Response Additional
    // ============================================================
    'LTCTE': {
        description: 'Thermal conductivity (electron).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Transport',
        wikiPage: 'LTCTE'
    },
    'LTETE': {
        description: 'Thermal transport (electron-electron).',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Transport',
        wikiPage: 'LTETE'
    },
    'LTHOMAS': {
        description: 'Thomas-Fermi screening.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Screening',
        wikiPage: 'LTHOMAS'
    },
    'LTDEP': {
        description: 'Temperature-dependent effective potential.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Phonons',
        wikiPage: 'LTDEP'
    },
    'LSYMGRAD': {
        description: 'Symmetrize gradient.',
        type: 'boolean',
        default: '.TRUE.',
        category: 'Symmetry',
        wikiPage: 'LSYMGRAD'
    },
    'LFOCKAEDFT': {
        description: 'Include Fock term in AE DFT.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Hybrid',
        wikiPage: 'LFOCKAEDFT'
    },
    'LVEL': {
        description: 'Velocity operator for optical properties.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Optics',
        wikiPage: 'LVEL'
    },
    'LDNEB': {
        description: 'Dynamical NEB.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NEB',
        wikiPage: 'LDNEB'
    },
    'LUSEW': {
        description: 'Use wavefunctions for transition state.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'NEB',
        wikiPage: 'LUSEW'
    },
    'LUSEVDW': {
        description: 'Use vdW-corrected energies.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Van der Waals',
        wikiPage: 'LUSEVDW'
    },
    'LCORR': {
        description: 'Harris-Foulkes correction to forces.',
        type: 'boolean',
        default: '.TRUE.',
        category: 'Basic',
        wikiPage: 'LCORR'
    },

    // ============================================================
    // Output Energy Terms
    // ============================================================
    'TOTEN': {
        description: 'Total free energy including entropy term (eV).',
        type: 'float',
        category: 'Output',
        wikiPage: 'TOTEN'
    },
    'E0': {
        description: 'Energy sigma->0 extrapolated (eV).',
        type: 'float',
        category: 'Output',
        wikiPage: 'TOTEN'
    },
    'EBANDS': {
        description: 'Band structure energy (eV).',
        type: 'float',
        category: 'Output',
        wikiPage: 'OUTCAR'
    },
    'EFERMI': {
        description: 'Fermi energy (eV).',
        type: 'float',
        category: 'Output',
        wikiPage: 'EFERMI'
    },

    // ============================================================
    // System Information
    // ============================================================
    'NKPTS': {
        description: 'Total number of k-points.',
        type: 'integer',
        category: 'K-points',
        wikiPage: 'KPOINTS'
    },
    'NIONS': {
        description: 'Total number of ions (atoms).',
        type: 'integer',
        category: 'System',
        wikiPage: 'POSCAR'
    },
    'NPLWV': {
        description: 'Total number of plane waves.',
        type: 'integer',
        category: 'System',
        wikiPage: 'OUTCAR'
    },
    'ALAT': {
        description: 'Lattice constant (Å).',
        type: 'float',
        category: 'System',
        wikiPage: 'POSCAR'
    },
    'A1': {
        description: 'First lattice vector.',
        type: 'float array',
        category: 'System',
        wikiPage: 'POSCAR'
    },
    'A2': {
        description: 'Second lattice vector.',
        type: 'float array',
        category: 'System',
        wikiPage: 'POSCAR'
    },
    'A3': {
        description: 'Third lattice vector.',
        type: 'float array',
        category: 'System',
        wikiPage: 'POSCAR'
    },
    'DIM': {
        description: 'Dimensionality of the system.',
        type: 'integer',
        category: 'System',
        wikiPage: 'OUTCAR'
    },
    'KPOINT': {
        description: 'Current k-point index.',
        type: 'integer',
        category: 'K-points',
        wikiPage: 'KPOINTS'
    },
    'KINTER': {
        description: 'K-point interpolation scheme.',
        type: 'integer',
        category: 'K-points',
        wikiPage: 'KINTER'
    },
    'NQ': {
        description: 'Number of q-points for response.',
        type: 'integer',
        category: 'Response',
        wikiPage: 'OUTCAR'
    },
    'QGRID': {
        description: 'Q-point grid specification.',
        type: 'integer array',
        category: 'Response',
        wikiPage: 'QGRID'
    },
    'IXMIN': {
        description: 'Minimum FFT index in x.',
        type: 'integer',
        category: 'Grid',
        wikiPage: 'OUTCAR'
    },
    'IXMAX': {
        description: 'Maximum FFT index in x.',
        type: 'integer',
        category: 'Grid',
        wikiPage: 'OUTCAR'
    },
    'SHIFTRED': {
        description: 'Shift reduction for band structures.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'K-points',
        wikiPage: 'SHIFTRED'
    },
    'EVENONLY': {
        description: 'Use only even k-points.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'K-points',
        wikiPage: 'EVENONLY'
    },
    'ODDONLY': {
        description: 'Use only odd k-points.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'K-points',
        wikiPage: 'ODDONLY'
    },
    'EVENONLYGW': {
        description: 'Use only even k-points for GW.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'EVENONLYGW'
    },
    'ODDONLYGW': {
        description: 'Use only odd k-points for GW.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'GW/BSE',
        wikiPage: 'ODDONLYGW'
    },
    'LMAXMP2': {
        description: 'Maximum L for MP2 calculations.',
        type: 'integer',
        category: 'GW/BSE',
        wikiPage: 'LMAXMP2'
    },

    // ============================================================
    // MD Additional Parameters
    // ============================================================
    'TAUPAR': {
        description: 'Parrinello-Rahman barostat time constant.',
        type: 'float',
        category: 'MD',
        wikiPage: 'TAUPAR'
    },
    'TEIN': {
        description: 'Initial kinetic energy for MD.',
        type: 'float',
        category: 'MD',
        wikiPage: 'TEIN'
    },
    'TELESCOPE': {
        description: 'Telescope sampling for response.',
        type: 'integer',
        category: 'Response',
        wikiPage: 'TELESCOPE'
    },
    'TURBO': {
        description: 'Turbo mode for faster calculations.',
        type: 'boolean',
        default: '.FALSE.',
        category: 'Algorithm',
        wikiPage: 'TURBO'
    },
    'VCA': {
        description: 'Virtual crystal approximation weights.',
        type: 'float array',
        category: 'VCA',
        wikiPage: 'VCA'
    },
    'RTIME': {
        description: 'Real-time TDDFT time step.',
        type: 'float',
        category: 'TDDFT',
        wikiPage: 'RTIME'
    },
    'MAXITERFT': {
        description: 'Maximum iterations for response function.',
        type: 'integer',
        category: 'Response',
        wikiPage: 'MAXITERFT'
    },
    'MCALPHA': {
        description: 'Alpha parameter for meta-GGA.',
        type: 'float',
        category: 'XC',
        wikiPage: 'MCALPHA'
    },
    'DEG_THRESHOLD': {
        description: 'Threshold for degeneracy detection.',
        type: 'float',
        category: 'Algorithm',
        wikiPage: 'DEG_THRESHOLD'
    },
};

// OUTCAR section patterns for navigation
const OUTCAR_SECTIONS = [
    { pattern: /^\s*Startparameter for this run/i, name: 'Startparameters' },
    { pattern: /^\s*Dimension of arrays/i, name: 'Dimension of arrays' },
    { pattern: /^\s*POSCAR:/i, name: 'POSCAR' },
    { pattern: /^\s*k-points\s+in\s+BZ/i, name: 'K-points' },
    { pattern: /^\s*Electronic Relaxation/i, name: 'Electronic Relaxation' },
    { pattern: /^\s*Ionic Relaxation/i, name: 'Ionic Relaxation' },
    { pattern: /^\s*DOS related values/i, name: 'DOS Parameters' },
    { pattern: /^\s*Iteration\s+(\d+)\s*\(/i, name: 'Iteration $1' },
    { pattern: /^\s*aborting loop/i, name: 'Loop Aborted' },
    { pattern: /^\s*FREE ENERGIE OF THE ION-ELECTRON SYSTEM/i, name: 'Free Energy' },
    { pattern: /^\s*VOLUME and BASIS/i, name: 'Volume and Basis' },
    { pattern: /^\s*TOTAL-FORCE/i, name: 'Total Force' },
    { pattern: /^\s*STRESS TENSOR/i, name: 'Stress Tensor' },
    { pattern: /^\s*EIGENVALUE/i, name: 'Eigenvalues' },
    { pattern: /^\s*E-fermi\s*:/i, name: 'Fermi Energy' },
    { pattern: /^\s*BRION:\s*g\(F\)/i, name: 'Ionic Convergence' },
    { pattern: /^\s*reached required accuracy/i, name: 'Converged!' },
    { pattern: /^\s*General timing and accounting/i, name: 'Timing Summary' },
    { pattern: /^\s*writing wavefunctions/i, name: 'Writing Wavefunctions' },
];

// Document symbol provider for OUTCAR navigation
class OutcarSymbolProvider implements vscode.DocumentSymbolProvider {
    public provideDocumentSymbols(
        document: vscode.TextDocument,
        _token: vscode.CancellationToken
    ): vscode.ProviderResult<vscode.DocumentSymbol[]> {
        const symbols: vscode.DocumentSymbol[] = [];
        const text = document.getText();
        const lines = text.split('\n');

        for (let i = 0; i < lines.length; i++) {
            const line = lines[i];

            for (const section of OUTCAR_SECTIONS) {
                const match = line.match(section.pattern);
                if (match) {
                    let name = section.name;
                    // Replace $1, $2, etc. with captured groups
                    if (match[1]) {
                        name = name.replace('$1', match[1]);
                    }

                    const range = new vscode.Range(i, 0, i, line.length);
                    const symbol = new vscode.DocumentSymbol(
                        name,
                        '',
                        vscode.SymbolKind.Function,
                        range,
                        range
                    );
                    symbols.push(symbol);
                    break; // Only match one pattern per line
                }
            }
        }

        return symbols;
    }
}

// Create hover provider for INCAR files
class VaspHoverProvider implements vscode.HoverProvider {
    public provideHover(
        document: vscode.TextDocument,
        position: vscode.Position,
        _token: vscode.CancellationToken
    ): vscode.ProviderResult<vscode.Hover> {
        const range = document.getWordRangeAtPosition(position, /[A-Za-z_][A-Za-z0-9_]*/);
        if (!range) {
            return null;
        }

        const word = document.getText(range).toUpperCase();
        const param = VASP_PARAMETERS[word];

        if (!param) {
            return null;
        }

        const markdown = new vscode.MarkdownString();
        markdown.isTrusted = true;

        // Parameter name and category
        markdown.appendMarkdown(`## ${word}\n\n`);
        markdown.appendMarkdown(`**Category:** ${param.category}\n\n`);

        // Description
        markdown.appendMarkdown(`${param.description}\n\n`);

        // Type and default
        markdown.appendMarkdown(`**Type:** \`${param.type}\`\n\n`);
        if (param.default) {
            markdown.appendMarkdown(`**Default:** \`${param.default}\`\n\n`);
        }

        // Possible values
        if (param.values && param.values.length > 0) {
            markdown.appendMarkdown(`**Values:** ${param.values.map(v => `\`${v}\``).join(', ')}\n\n`);
        }

        // Link to VASP wiki
        const wikiUrl = `${VASP_WIKI_BASE}${param.wikiPage}`;
        markdown.appendMarkdown(`---\n\n`);
        markdown.appendMarkdown(`[Open VASP Wiki documentation](${wikiUrl})\n`);

        return new vscode.Hover(markdown, range);
    }
}

// Provide completions for INCAR parameters
class VaspCompletionProvider implements vscode.CompletionItemProvider {
    public provideCompletionItems(
        document: vscode.TextDocument,
        position: vscode.Position,
        _token: vscode.CancellationToken,
        _context: vscode.CompletionContext
    ): vscode.ProviderResult<vscode.CompletionItem[]> {
        const linePrefix = document.lineAt(position).text.substring(0, position.character);

        // Only complete at the start of a line or after whitespace
        if (!/^[\s]*$/.test(linePrefix) && !/[\s=]$/.test(linePrefix)) {
            // Check if we're typing a parameter name
            if (!/^[\s]*[A-Za-z_][A-Za-z0-9_]*$/.test(linePrefix)) {
                return [];
            }
        }

        const completions: vscode.CompletionItem[] = [];

        for (const [name, param] of Object.entries(VASP_PARAMETERS)) {
            const item = new vscode.CompletionItem(name, vscode.CompletionItemKind.Property);
            item.detail = `[${param.category}] ${param.type}`;
            item.documentation = new vscode.MarkdownString(
                `${param.description}\n\n` +
                (param.default ? `**Default:** \`${param.default}\`\n\n` : '') +
                `[VASP Wiki](${VASP_WIKI_BASE}${param.wikiPage})`
            );

            // Add snippet for boolean parameters
            if (param.type === 'boolean') {
                item.insertText = new vscode.SnippetString(`${name} = \${1|.TRUE.,.FALSE.|}`);
            } else if (param.values && param.values.length > 0) {
                const choices = param.values.join(',');
                item.insertText = new vscode.SnippetString(`${name} = \${1|${choices}|}`);
            } else {
                item.insertText = new vscode.SnippetString(`${name} = $1`);
            }

            completions.push(item);
        }

        return completions;
    }
}

export function activate(context: vscode.ExtensionContext) {
    console.log('VASP Syntax extension is now active');

    // Register hover provider for INCAR and OUTCAR files
    const hoverProvider = new VaspHoverProvider();
    context.subscriptions.push(
        vscode.languages.registerHoverProvider('incar', hoverProvider),
        vscode.languages.registerHoverProvider('outcar', hoverProvider)
    );

    // Register document symbol provider for OUTCAR navigation
    const outcarSymbolProvider = new OutcarSymbolProvider();
    context.subscriptions.push(
        vscode.languages.registerDocumentSymbolProvider('outcar', outcarSymbolProvider)
    );

    // Register completion provider for INCAR files
    const completionProvider = new VaspCompletionProvider();
    context.subscriptions.push(
        vscode.languages.registerCompletionItemProvider('incar', completionProvider)
    );

    // Register a command to open VASP wiki for the word under cursor
    context.subscriptions.push(
        vscode.commands.registerCommand('vasp.openWiki', () => {
            const editor = vscode.window.activeTextEditor;
            if (!editor) {
                return;
            }

            const position = editor.selection.active;
            const range = editor.document.getWordRangeAtPosition(position, /[A-Za-z_][A-Za-z0-9_]*/);
            if (!range) {
                return;
            }

            const word = editor.document.getText(range).toUpperCase();
            const param = VASP_PARAMETERS[word];

            if (param) {
                const wikiUrl = `${VASP_WIKI_BASE}${param.wikiPage}`;
                vscode.env.openExternal(vscode.Uri.parse(wikiUrl));
            } else {
                vscode.window.showInformationMessage(`No VASP documentation found for '${word}'`);
            }
        })
    );
}

export function deactivate() {}
