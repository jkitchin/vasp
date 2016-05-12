import glob
import os


def task_test():
    return {'actions': ["nosetests"]}


# These seem to stop on the first error.
def task_flakes():
    """flakes on targets"""
    for py in glob.glob("vasp/*.py"):
        yield {'basename': 'Flakes',
               'name': py,
               'actions': ["ls {}".format(py)],
               'file_dep': [py],
               'verbosity': 2}


def task_pep8():
    for py in glob.glob("vasp/*.py"):
        yield {'basename': 'pep8',
               'name': py,
               'actions': [['pep8',  py]],
               'file_dep': [py],
               'verbosity': 2}
