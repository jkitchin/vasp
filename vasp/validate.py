"""Validation functions for Vasp keywords.

Each function should have the signature func(calc, val) and it should use
exceptions or assertions to validate. Each function should have a brief
docstring. The first line will be used as a tooltip in Emacs. A command will
give access to the full docstring. It is encouraged to put URLs to full
documentation, as they will be clickable in Emacs.

"""
import types
from vasp import Vasp


def ediffg(calc, val):
    """EDIFFG defines the break condition for the ionic relaxation loop.

    If EDIFFG < 0, it defines a force criteria.
    """
    assert isinstance(val, float) or val == 0


def encut(calc, val):
    """Planewave cutoff in eV."""
    assert val > 0, 'encut must be greater than zero.'
    assert isinstance(val, int) or isinstance(val, float),\
        'encut should be an int or float'


def ibrion(calc, val):
    """IBRION determines how the ions are updated and moved."""
    assert val in [-1, 0, 1, 2, 3, 5, 6, 7, 8, 44]


def isif(calc, val):
    """ISIF determines which DOF (ions, cell volume, cell shape) are allowed to
    change.

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
    """ISMEAR determines how the partial occupancies $ f_{n{\bf k}}$ are set for
    each orbital.

    """
    assert val in [-5, -4, -3, -2, 0, 1, 2]


def ispin(calc, val):
    """Control spin-polarization.

    1 - default, no spin polarization
    2 - spin-polarization.

    """
    assert val in [1, 2], "ispin should be 1 or 2"
    if val == 2:
        assert 'magmom' in calc.parameters, "magmom is not set."
        assert len(calc.parameters['magmom']) == len(calc.get_atoms()),\
                   "len(magmom) != len(atoms)"


def nsw(calc, val):
    """NSW sets the maximum number of ionic steps."""
    assert isinstance(val, int)


def prec(calc, val):
    """Precision."""
    assert val in ['Low', 'Medium', 'High', 'Normal', 'Accurate', 'Single']


def sigma(calc, val):
    """SIGMA determines the width of the smearing in eV."""
    assert isinstance(val, float)
    assert val > 0


def xc(calc, val):
    """Set exchange-correlation functional."""
    assert val in Vasp.xc_defaults.keys(), \
        "xc not in {}.".format(Vasp.xc_defaults.keys())


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
