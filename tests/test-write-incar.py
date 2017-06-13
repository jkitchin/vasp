from nose import *
from vasp import *
import os
from ase import Atom, Atoms

def setup_func():
    "set up test fixtures"
    if os.path.exists('INCAR'):
        os.unlink('INCAR')


def teardown_func():
    "tear down test fixtures"
    if os.path.exists('INCAR'):
        os.unlink('INCAR')


@with_setup(setup_func, teardown_func)
def test0():
    "check ispin"

    atoms = Atoms([Atom('O', [4, 5, 5], magmom=1),
                   Atom('C', [5, 5, 5], magmom=2),
                   Atom('O', [6, 5, 5], magmom=3)],
                   cell=(10, 10, 10))

    calc = Vasp('vasp',
                ispin=2,
                atoms=atoms)

    calc.write_incar('INCAR')

    incar = calc.read_incar('INCAR')

    assert incar['lorbit'] == 11
    assert incar['ispin'] == 2
    assert incar['magmom'] == [1.0, 3.0, 2.0]

    # Now delete it.
    calc.set(ispin=None)
    calc.write_incar('INCAR')

    incar = calc.read_incar('INCAR')
    assert 'ispin' not in incar
    assert 'lorbit' not in incar
    assert 'magmom' not in incar


@with_setup(setup_func, teardown_func)
def test1():
    atoms = Atoms([Atom('O', [4, 5, 5], magmom=1),
                   Atom('C', [5, 5, 5], magmom=2),
                   Atom('O', [6, 5, 5], magmom=3)],
                   cell=(10, 10, 10))

    calc = Vasp('vasp',
                ispin=2,
                nsw=1,
                atoms=atoms)

    calc.write_incar('INCAR')

    incar = calc.read_incar('INCAR')
    print(incar)
    assert incar['lwave'] is False
    assert incar['nsw'] == 1

    calc.set(nsw=0)
    calc.write_incar('INCAR')
    incar = calc.read_incar('INCAR')
    assert incar['lwave'] is False
    assert incar['nsw'] == 0

    calc.set(nsw=None)
    calc.write_incar('INCAR')
    incar = calc.read_incar('INCAR')
    assert incar['lwave'] is False
    assert 'nsw' not in incar


# default settings
@with_setup(setup_func, teardown_func)
def test2():
    "check default sigma"

    atoms = Atoms([Atom('O', [4, 5, 5], magmom=1),
                   Atom('C', [5, 5, 5], magmom=2),
                   Atom('O', [6, 5, 5], magmom=3)],
                   cell=(10, 10, 10))

    calc = Vasp('vasp',
                atoms=atoms)

    calc.write_incar('INCAR')

    incar = calc.read_incar('INCAR')

    assert incar['sigma'] == 0.1


@with_setup(setup_func, teardown_func)
def test2():
    "check default sigma"

    atoms = Atoms([Atom('O', [4, 5, 5], magmom=1),
                   Atom('C', [5, 5, 5], magmom=2),
                   Atom('O', [6, 5, 5], magmom=3)],
                   cell=(10, 10, 10))

    calc = Vasp('vasp',
                sigma=0.01,
                atoms=atoms)

    calc.write_incar('INCAR')

    incar = calc.read_incar('INCAR')

    assert incar['sigma'] == 0.01


@with_setup(setup_func, teardown_func)
def test3():
    "check rwigs"

    atoms = Atoms([Atom('O', [4, 5, 5], magmom=1),
                   Atom('C', [5, 5, 5], magmom=2),
                   Atom('O', [6, 5, 5], magmom=3)],
                   cell=(10, 10, 10))

    calc = Vasp('vasp',
                rwigs={'C': 1.5, 'O': 1.0},
                atoms=atoms)

    calc.write_incar('INCAR')

    incar = calc.read_incar('INCAR')
    assert incar['rwigs'] == [1.0, 1.5]
