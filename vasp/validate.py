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


def ibrion(calc, val):
    """IBRION determines the algorithm to update geometry during relaxtion. (int)

    http://cms.mpi.univie.ac.at/wiki/index.php/IBRION

    """
    assert val in [-1, 0, 1, 2, 3, 5, 6, 7, 8, 44]


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


def kpts(calc, val):
    """Sets k-points. Not a Vasp keyword. (list or tuple)"""
    assert isinstance(val, list) or isinstance(val, tuple)


def kpts_nintersections(calc, val):
    """Triggers line mode in KPOINTS for bandstructure calculations. (int)

    http://cms.mpi.univie.ac.at/vasp/vasp/Strings_k_points_bandstructure_calculations.html
    """
    assert isinstance(val, int)


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
    assert val > calc.get_valence_electrons() / 2, \
        'nbands = {} which is less than {}.'.format(val, calc.get_valence_electrons() / 2)


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
