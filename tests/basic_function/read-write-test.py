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
                    'atom-function-test')


def teardown_func():
    # Remove test directory after finishing

    if os.path.exists('atom-function-test'):
        shutil.rmtree('atom-function-test')


@with_setup(setup_func, teardown_func)
def test():
    # This calculation is stored in "premade_calculations/co2"
    atoms = Atoms([Atom('O', [4, 5, 5]),
                   Atom('C', [5, 5, 5]),
                   Atom('O', [6, 5, 5])],
                  cell=(10, 10, 10))

    # Initiate calculator with atoms object
    calc = Vasp('atom-function-test', atoms=atoms)

    # Call various properties
    atoms = calc.get_atoms()
    atoms.get_potential_energy()
    atoms.get_cell()
    atoms.get_forces()
    atoms.get_stress()

    # Initiate without atoms object
    calc = Vasp('atom-function-test')
    atoms = calc.get_atoms()
    atoms.get_potential_energy()
    atoms.get_cell()
    atoms.get_forces()
    atoms.get_stress()


if __name__ == '__main__':
    test()
