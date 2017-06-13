from nose import with_setup
from vasp import Vasp
import os
from ase import Atom, Atoms
from ase.db import connect


def setup_func():
    "set up test fixtures"
    if os.path.exists('DB.db'):
        os.unlink('DB.db')


def teardown_func():
    "tear down test fixtures"
    if os.path.exists('DB.db'):
        os.unlink('DB.db')


@with_setup(setup_func, teardown_func)
def test0():
    atoms = Atoms([Atom('O', [4, 5, 5], magmom=1),
                   Atom('C', [5, 5, 5], magmom=2),
                   Atom('O', [6, 5, 5], magmom=3)],
                   cell=(10, 10, 10))

    calc = Vasp('vasp',
                sigma=0.01,
                atoms=atoms)

    calc.write_db(fname='DB.db', append=False)

    with connect('DB.db') as con:
        data = con.get(1).data
        assert data['parameters']['sigma'] == 0.01

    print('done')
