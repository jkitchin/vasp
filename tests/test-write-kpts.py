from nose import *
from vasp import *
import os


def setup_func():
    "set up test fixtures"
    if os.path.exists('KPOINTS'):
        os.unlink('KPOINTS')


def teardown_func():
    "tear down test fixtures"
    if os.path.exists('KPOINTS'):
        os.unlink('KPOINTS')


@with_setup(setup_func, teardown_func)
def test0():
    "write automatic with Monkhorstpack grid"
    calc = Vasp('vasp', kpts=[4, 4, 4])
    calc.write_kpoints('KPOINTS')
    result = open('KPOINTS', 'r').read()
    assert result == '''\
KPOINTS created by Atomic Simulation Environment
0
Monkhorst-Pack
4 4 4
0.0 0.0 0.0
'''


@with_setup(setup_func, teardown_func)
def test1():
    "write automatic with Monkhorstpack grid"
    calc = Vasp('vasp',
                kpts_nintersections=10,
                reciprocal=True,
                kpts=[[0.5, 0.5, 0.0],
                      [0, 0, 0],
                      [0, 0, 0],
                      [0.5, 0.5, 0.5]])
    calc.write_kpoints('KPOINTS')
    result = open('KPOINTS', 'r').read()
    assert result == '''\
KPOINTS created by Atomic Simulation Environment
10
Line-mode
Reciprocal
0.5 0.5 0.0
0 0 0
0 0 0
0.5 0.5 0.5
'''
