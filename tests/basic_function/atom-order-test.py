#!/usr/bin/env python
from nose.tools import with_setup
from ase import Atoms, Atom
from vasp import Vasp
import shutil
import os


def setup_func():
    # Calculations can take awhile, so we
    # simulate a finished calc by copying it in.

    shutil.copytree('premade_calculations/co2',
                    './atom-order-test')


def teardown_func():
    # Remove test directory after finishing

    if os.path.exists('atom-order-test'):
        shutil.rmtree('atom-order-test')


@with_setup(setup_func, teardown_func)
def test():
    # This calculation is stored in "premade_calculations/co2"
    atoms0 = Atoms([Atom('O', [4, 5, 5]),
                    Atom('C', [5, 5, 5]),
                    Atom('O', [6, 5, 5])],
                   cell=(10, 10, 10))

    # Test that atoms object agrees with calculator
    calc = Vasp('atom-order-test',
                xc='PBE',
                encut=350,
                ediff=0.1,
                nsw=0,
                atoms=atoms0)

    # Standard call
    atoms = calc.get_atoms()
    atoms.get_potential_energy()

    # Check the atom order
    if not atoms.get_chemical_symbols() == \
       atoms0.get_chemical_symbols():
        raise Exception('Atom order not conserved')

    # Now call the calculator without assigning atoms
    calc2 = Vasp('atom-order-test')
    atoms = calc2.get_atoms()

    # Check the atom order
    if not atoms.get_chemical_symbols() == \
       atoms0.get_chemical_symbols():
        raise Exception('Atom order not conserved')


if __name__ == '__main__':
    test()
