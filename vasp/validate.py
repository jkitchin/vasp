"""Validation functions for Vasp keywords.

Each function should have the signature func(calc, val) and it should
use exceptions or assertions to validate. Each function should have a
brief docstring. The first line will be used as a tooltip in Emacs. An
Emacs command will give access to the full docstring. It is encouraged
to put URLs to full documentation, as they will be clickable in Emacs.

http://cms.mpi.univie.ac.at/wiki/index.php/Category:INCAR

"""
import types
import ase
from ase.utils import basestring
import warnings


def algo(calc, val):
    """ specify the electronic minimisation algorithm (as of VASP.4.5) (string)

    http://cms.mpi.univie.ac.at/wiki/index.php/ALGO

    """
    assert isinstance(val, str)
    assert val in ["Normal", "VeryFast", "Fast" , "Conjugate", "All", "Damped",
                   "Subrot", "Eigenval", "None", "Nothing", "CHI", "GW0",
                   "GW", "scGW0", "scGW"]


def atoms(calc, val):
    """The Atoms object. (ase.atoms.Atoms or a list of them for an NEB)."""
    assert isinstance(val, ase.atoms.Atoms) or isinstance(val, list)


def ediff(calc, val):
    """EDIFF specifies the global break condition for the electronic loop. (float)

    http://cms.mpi.univie.ac.at/wiki/index.php/EDIFF
    """
    assert isinstance(val, float) or val == 0


def ediffg(calc, val):
    """EDIFFG defines the break condition for the ionic relaxation loop. (float)

    If EDIFFG < 0, it defines a force criteria.

    http://cms.mpi.univie.ac.at/wiki/index.php/EDIFFG
    """
    assert isinstance(val, float) or val == 0


def encut(calc, val):
    """Planewave cutoff in eV. (float)

    http://cms.mpi.univie.ac.at/wiki/index.php/ENCUT
    """
    assert val > 0, 'encut must be greater than zero.'
    assert (isinstance(val, int)
            or isinstance(val, long)
            or isinstance(val, float)),\
        ('encut should be an int or float.'
         ' You provided {} ({}).'.format(val, type(val)))


def ialgo(calc, val):
    """IALGO selects the algorithm used to optimize the orbitals. (int)

    http://cms.mpi.univie.ac.at/wiki/index.php/IALGO

    """
    print('You are advised to use the algo key instead of ialgo.')
    assert isinstance(val, int)
    assert val in [-1, 2, 3, 4, 5, 6, 7, 8,
                   15, 16, 17, 18,
                   28,
                   38,
                   44, 45, 46, 47, 48,
                   53, 54, 55, 56, 57, 58]


def ibrion(calc, val):
    """IBRION determines the algorithm to update geometry during relaxtion. (int)

    http://cms.mpi.univie.ac.at/wiki/index.php/IBRION

    """
    assert val in [-1, 0, 1, 2, 3, 5, 6, 7, 8, 44]


def icharg(calc, val):
    """Determines how VASP constructs the initial charge density. (int)

    0 - calculate from initial wave functions
    1 - read from the CHGCAR
    2 - (default) Superposition of atomic charge densities
    11 - for band-structure plots

    http://cms.mpi.univie.ac.at/wiki/index.php/ICHARG

    """
    assert isinstance(val, int)


def images(calc, val):
    """The number of images not counting the end-points for an NEB. (int)

    http://cms.mpi.univie.ac.at/wiki/index.php/IMAGES

    """
    assert isinstance(val, int)
    assert val == len(calc.neb) - 2


def isif(calc, val):
    """ISIF determines what is changed during relaxations. (int)

    | ISIF | calculate | calculate        | relax | change     | change      |
    |      | force     | stress tensor    | ions  | cell shape | cell volume |
    |------+-----------+------------------+-------+------------+-------------|
    |    0 | yes       | no               | yes   | no         | no          |
    |    1 | yes       | trace only $ ^*$ | yes   | no         | no          |
    |    2 | yes       | yes              | yes   | no         | no          |
    |    3 | yes       | yes              | yes   | yes        | yes         |
    |    4 | yes       | yes              | yes   | yes        | no          |
    |    5 | yes       | yes              | no    | yes        | no          |
    |    6 | yes       | yes              | no    | yes        | yes         |
    |    7 | yes       | yes              | no    | no         | yes         |

    """
    assert val in [0, 1, 2, 3, 4, 5, 6, 7]


def ismear(calc, val):
    """ISMEAR determines how the partial occupancies are set (int).

    http://cms.mpi.univie.ac.at/wiki/index.php/ISMEAR

    """
    assert val in [-5, -4, -3, -2, -1, 0, 1, 2]


def ispin(calc, val):
    """Control spin-polarization. (int)

    1 - default, no spin polarization
    2 - spin-polarization.

    http://cms.mpi.univie.ac.at/wiki/index.php/ISPIN

    """
    assert val in [1, 2], "ispin should be 1 or 2"
    if val == 2:
        assert 'magmom' in calc.parameters, "magmom is not set."
        assert len(calc.parameters['magmom']) == len(calc.get_atoms()),\
                   "len(magmom) != len(atoms)"


def isym(calc, val):
    """ISYM determines the way VASP treats symmetry.

    http://cms.mpi.univie.ac.at/wiki/index.php/ISYM
    """
    assert val in [-1, 0, 1, 2, 3]


def ldau(calc, val):
    """ LDAU switches on the L(S)DA+U. (bool)

    http://cms.mpi.univie.ac.at/wiki/index.php/LDAU
    """
    assert val in [True, False, None]


def ldau_luj(calc, val):
    """Dictionary of DFT+U parameters.

    ldau_luj={'Mn': {'L': 2, 'U': 0.0, 'J': 0.0},
              'O': {'L': -1, 'U': 0.0, 'J': 0.0}},
    """
    assert isinstance(val, dict)
    # this may not be the case for site-specific U. I think we need
    # setups for that.
    assert len(val.keys()) == len(set([a.symbol for a in calc.get_atoms()]))


def ldauprint(calc, val):
    """LDAUPRINT controls the verbosity of the L(S)DA+U routines. (int)

    LDAUPRINT=0: silent.
    LDAUPRINT=1: Write occupancy matrix to the OUTCAR file.
    LDAUPRINT=2: same as LDAUPRINT=1, plus potential matrix dumped to stdout.

    http://cms.mpi.univie.ac.at/wiki/index.php/LDAUPRINT
    """
    assert val in [0, 1, 2]


def ldautype(calc, val):
    """LDAUTYPE specifies which type of L(S)DA+U approach will be used. (int)

    LDAUTYPE=1: The rotationally invariant LSDA+U introduced by
    Liechtenstein et al.[1]

    LDAUTYPE=2: The simplified (rotationally invariant) approach to
    the LSDA+U, introduced by Dudarev et al.[2]

    1. A. I. Liechtenstein, V. I. Anisimov and J. Zaane, Phys. Rev. B
    52, R5467 (1995).

    2. S. L. Dudarev, G. A. Botton, S. Y. Savrasov, C. J. Humphreys
    and A. P. Sutton, Phys. Rev. B 57, 1505 (1998).

    http://cms.mpi.univie.ac.at/wiki/index.php/LDAUTYPE

    """


def lmaxmix(calc, val):
    """LMAXMIX the max l-quantum number the charge densities used. (int)

    Mostly used for DFT+U.
    4 for d-electrons (or 6 for f-elements)

    http://cms.mpi.univie.ac.at/wiki/index.php/LMAXMIX
    """
    assert val in [2, 4, 6]


def kpts(calc, val):
    """Sets k-points. Not a Vasp keyword. (list)"""
    assert isinstance(val, list), 'kpts should be a list.'


def kpts_nintersections(calc, val):
    """Triggers line mode in KPOINTS for bandstructure calculations. (int)

    http://cms.mpi.univie.ac.at/vasp/vasp/Strings_k_points_bandstructure_calculations.html
    """
    assert isinstance(val, int)


def kspacing(calc, val):
    """KSPACING determines the number of k-points if the KPOINTS file is not present (float).

    http://cms.mpi.univie.ac.at/vasp/vasp/KSPACING_tag_KGAMMA_tag.html
    """
    assert(isinstance(val, float))


def lcharg(calc, val):
    """LCHARG determines whether CHGCAR and CHG are written. (boolean)

    http://cms.mpi.univie.ac.at/wiki/index.php/LCHARG

    """
    assert val in [True, False]


def lorbit(calc, val):
    """Determines whether the PROCAR or PROOUT files are written.

    http://cms.mpi.univie.ac.at/wiki/index.php/LORBIT

    """
    if val < 10:
        assert 'rwigs' in calc.parameters
    assert isinstance(val, int)


def lreal(calc, val):
    """LREAL determines whether the projection operators are evaluated in real-space or in reciprocal space. (boolean)

    http://cms.mpi.univie.ac.at/wiki/index.php/LREAL

    """
    assert val in [True, False, 'On', 'Auto', 'O', 'A']


def lwave(calc, val):
    """LWAVE determines whether the WAVECAR is written. (Boolean)

    http://cms.mpi.univie.ac.at/wiki/index.php/LWAVE

    """
    assert val in [True, False]


def magmom(calc, val):
    """MAGMOM Specifies the initial magnetic moment for each atom. (list)

    http://cms.mpi.univie.ac.at/wiki/index.php/MAGMOM

    """
    assert isinstance(val, list),\
        'Got {} for magmom. Should be a list.'.format(val)

    assert len(val) == len(calc.atoms)


def nbands(calc, val):
    """NBANDS determines the actual number of bands in the calculation. (int)

    http://cms.mpi.univie.ac.at/wiki/index.php/NBANDS

    """
    assert isinstance(val, int) or isinstance(val, long)

    s = 'nbands = {} which is less than {}.'
    assert val > calc.get_valence_electrons() / 2, \
        s.format(val, calc.get_valence_electrons() / 2)


def nelm(calc, val):
    """NELM sets the maximum number of electronic SC (selfconsistency) steps which may be performed. (int)

    http://cms.mpi.univie.ac.at/wiki/index.php/NELM
    """

    assert isinstance(val, int)


def nupdown(calc, val):
    """NUPDOWN = difference between number of spin up and down electrons.

    This fixes the bulk magnetic moment.
    The VASP manual specifies this should be an integer, but it
    appears floats work too.

    http://cms.mpi.univie.ac.at/vasp/vasp/NUPDOWN.html

    """
    assert isinstance(val, int) or isinstance(val, float)


def nsw(calc, val):
    """NSW sets the maximum number of ionic steps. (int)

    http://cms.mpi.univie.ac.at/wiki/index.php/NSW

    """
    assert isinstance(val, int)


def pp(calc, val):
    """Determines where POTCARS are retrieved from. (string)"""
    assert val in ['PBE', 'LDA', 'GGA']


def prec(calc, val):
    """Specifies the Precision-mode. (string)

    http://cms.mpi.univie.ac.at/wiki/index.php/PREC
    """
    assert val in ['Low', 'Medium', 'High', 'Normal', 'Accurate', 'Single']


def reciprocal(calc, val):
    """Specifies reciprocal coordinates in KPOINTS. (boolean)

    Not a Vasp keyword."""
    assert val in [True, False]


def rwigs(calc, val):
    """RWIGS specifies the Wigner-Seitz radius for each atom type. (list)

    in vasp.py you enter a dictionary of {sym: radius}.

    http://cms.mpi.univie.ac.at/wiki/index.php/RWIGS
    """
    assert isinstance(val, list)
    assert calc.parameters.get('lorbit', 0) < 10, \
        'lorbit >= 10, rwigs is ignored.'


def setups(calc, val):
    """Sets up special setups for the POTCARS (list of (symbol/int, suffix)).

    The first element of each pair of the list is either an integer
    index of the atom for the special setup, or a chemical symbol for
    all atoms of that type. The second element of the pair is a suffix
    to be appended to the symbol. For example, to use the O_s potcar
    set setups to: [['O', '_s']].

    This is not a vasp keyword.

    """
    assert isinstance(val, list)
    for s, suffix in val:
        assert isinstance(s, int) or isinstance(s, basestring)
        assert isinstance(suffix, basestring)


def sigma(calc, val):
    """SIGMA determines the width of the smearing in eV. (float)"""
    assert isinstance(val, float)
    assert val > 0


def spring(calc, val):
    """The Spring constant  in the elastic band method. -5 = NEB.

    http://cms.mpi.univie.ac.at/wiki/index.php/SPRING
    """
    assert isinstance(val, int) or isinstance(val, float)
    if calc.parameters.get('ibrion') not in [1, 3]:
        warnings.warn('ibrion should be 1 or 3.')


def xc(calc, val):
    """Set exchange-correlation functional. (string)"""
    import vasp
    assert val.lower() in vasp.Vasp.xc_defaults.keys(), \
        "xc ({}) not in {}.".format(val, vasp.Vasp.xc_defaults.keys())


def keywords():
    """Return list of keywords we validate.

    Returns a lisp list for Emacs.

    """
    import validate

    f = [validate.__dict__.get(a) for a in dir(validate)
         if isinstance(validate.__dict__.get(a), types.FunctionType)]

    names = [x.__name__ for x in f]
    names.remove('keywords')

    return "(" + ' '.join(['"{}"'.format(x) for x in names]) + ")"


def keyword_alist():
    """Returns an alist of (keyword . "first doc string").

    Returns the alist for use in Emacs.

    """
    import validate
    f = [validate.__dict__.get(a) for a in dir(validate)
         if isinstance(validate.__dict__.get(a), types.FunctionType)]

    names = [x.__name__ for x in f]
    names.remove('keywords')
    names.remove('keyword_alist')

    docstrings = [validate.__dict__[name].__doc__.split('\n')[0]
                  for name in names]

    cons_cells = ["(\"{}\"  \"{}\")".format(key, doc)
                  for key, doc in zip(names, docstrings)]
    return "(" + "".join(cons_cells) + ")"
