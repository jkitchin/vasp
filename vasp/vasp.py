"""ASE calculator for Vasp."""

# https://wiki.fysik.dtu.dk/ase/_modules/ase/calculators/calculator.html#Calculator

# https://gitlab.com/ase/ase/blob/master/ase/calculators/abinit.py

# ###################################################################
# Logger for handling information, warning and debugging
# ###################################################################
import logging
import sys

log = logging.getLogger('Vasp')
log.setLevel(logging.CRITICAL)
handler = logging.StreamHandler()
if sys.version_info < (2, 5):  # no funcName in python 2.4
    formatstring = ('%(levelname)-10s '
                    'lineno: %(lineno)-4d %(message)s')
else:
    formatstring = ('%(levelname)-10s function: %(funcName)s '
                    'lineno: %(lineno)-4d %(message)s')
formatter = logging.Formatter(formatstring)
handler.setFormatter(formatter)
log.addHandler(handler)


# The core class is defined here
from vasp_core import Vasp

# These modules monkeypatch the Vasp class
import writers
import readers
import setters
import runner
