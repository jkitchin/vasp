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

import sys

# I want all the functions wrapped in a try block so we can use our
# own exception handler.
def tryit(func):
    """Decorator to wrap function in try block with Vasp exception handler."""
    def inner(self, *args, **kwargs):
        try:
            return func(self, *args, **kwargs)
        except Exception:
            if self.exception_handler is not None:
                return self.exception_handler(*sys.exc_info())
            else:
                raise

    inner.__name__ = func.__name__
    inner.__doc__ = func.__doc__
    return inner

from ase.calculators.calculator import Calculator,\
    FileIOCalculator


for attr in Vasp.__dict__:
    if callable(getattr(Vasp, attr)):
        setattr(Vasp, attr, tryit(getattr(Vasp, attr)))

for attr in Calculator.__dict__:
    if callable(getattr(Calculator, attr)):
        setattr(Calculator, attr, tryit(getattr(Calculator, attr)))

for attr in FileIOCalculator.__dict__:
    if callable(getattr(FileIOCalculator, attr)):
        setattr(FileIOCalculator, attr, tryit(getattr(FileIOCalculator, attr)))
