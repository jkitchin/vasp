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
import getters
import setters
import vib
import neb
import serialize
import runner
import bader


def tryit(func):
    """Decorator to wrap function in try block with Vasp exception handler."""
    def inner(self, *args, **kwargs):
        # we have to check for this attr because sometimes the
        # calculator is a singlepoint calculator which doesn't have
        # it. and also for when we write to the ase-db
        if (hasattr(self, 'debug') and self.debug is not None
            or (not hasattr(self, 'exception_handler'))):
            return func(self, *args, **kwargs)
        else:
            try:
                return func(self, *args, **kwargs)
            except Exception:
                if self.exception_handler is not None:
                    return self.exception_handler(self, *sys.exc_info())
                else:
                    raise

    inner.__name__ = func.__name__
    if hasattr(func, 'im_func'):
        template = '''

        Wrapped in vasp.tryit.

        Defined at [[{0}::{1}]].'''.format(func.im_func.func_code.co_filename,
                                           func.im_func.func_code.co_firstlineno)
    else:
        template = '''

        Wrapped in vasp.tryit.'''

    if func.__doc__ is not None:
        inner.__doc__ = func.__doc__ + template
    else:
        inner.__doc__ = 'Wrapped in vasp.tryit.'
    return inner

from ase.calculators.calculator import Calculator,\
    FileIOCalculator

import inspect

# We avoid decorating class methods. It seems to break them.
for attr in Vasp.__dict__:
    f = getattr(Vasp, attr)
    if inspect.ismethod(f) and f.__self__ is not Vasp:
        setattr(Vasp, attr, tryit(getattr(Vasp, attr)))

for attr in Calculator.__dict__:
    f = getattr(Vasp, attr)
    if inspect.ismethod(f) and f.__self__ is not Vasp:
        setattr(Calculator, attr, tryit(getattr(Calculator, attr)))

for attr in FileIOCalculator.__dict__:
    f = getattr(Vasp, attr)
    if inspect.ismethod(f) and f.__self__ is not Vasp:
        setattr(FileIOCalculator, attr,
                tryit(getattr(FileIOCalculator, attr)))
