"""Writer functions for vasp.py

Functions that write files: INCAR, POSCAR, POTCAR, KPOINTS

These are separated out by design to keep vasp.py small. Each function is
monkey-patched onto the Vasp class as if it were defined in vasp.py.

"""
import os
import numpy as np
import vasp
from monkeypatch import monkeypatch_class
from ase.calculators.calculator import FileIOCalculator


@monkeypatch_class(vasp.Vasp)
def write_input(self, atoms=None, properties=None, system_changes=None):
    """Writes all input files required for a calculation."""
    # this creates the directory if needed
    FileIOCalculator.write_input(self, atoms, properties, system_changes)

    if 'spring' not in self.parameters:  # do not write if NEB
        self.write_poscar()
    self.write_incar()
    self.write_kpoints()
    self.write_potcar()
    self.write_db()


@monkeypatch_class(vasp.Vasp)
def write_db(self, atoms=None, fname=None, data=None, **kwargs):
    """Write the DB file.

    atoms can be any atoms object, defaults to self.get_atoms().
    fname can be anything, defaults to self.directory/DB.db

    data is a dictionary of data to store.

    kwargs is key=value pairs to store with the atoms.

    Existing data and kwargs are preserved. You can delete kwargs by
    setting them to None. You can delete data by setting the key to
    None in the data dictionary.

    Only row 1 should be in this database.

    """
    from ase.db import connect

    if fname is None:
        fname = os.path.join(self.directory, 'DB.db')

    fdata = {'resort': self.resort,
             'parameters': self.parameters,
             'ppp_list': self.ppp_list}
    fkv = {'path': self.directory}

    # get current data and keywords
    if os.path.exists(fname):

        with connect(fname) as con:
            try:
                temp_atoms = con.get_atoms(id=1,
                                           add_additional_information=True)
                fdata.update(temp_atoms.info['data'])
                fkv.update(temp_atoms.info['key_value_pairs'])
            except KeyError:
                pass

    # Update fdata from input args. None removes keywords and data
    # elements
    if data is not None:
        for key, val in data.iteritems():
            if val is None and key in fdata:
                del fdata[key]
            else:
                fdata.update({key: val})

    # update key-value pairs from input args
    for key, val in kwargs.iteritems():
        # we use None to delete keys
        if val is None and key in fkv:
            del fkv[key]
        elif val is not None:
            fkv.update({key: val})

    if atoms is None:
        atoms = self.get_atoms()

    # TODO: NEB? should contain images?
    # write out in non-append mode.
    with connect(fname, append=False) as con:
        con.write(atoms,
                  data=fdata,
                  **fkv)


@monkeypatch_class(vasp.Vasp)
def write_poscar(self, fname=None):
    """Write the POSCAR file."""
    if fname is None:
        fname = os.path.join(self.directory, 'POSCAR')

    from ase.io.vasp import write_vasp
    write_vasp(fname,
               self.atoms_sorted,
               symbol_count=self.symbol_count)


@monkeypatch_class(vasp.Vasp)
def write_incar(self, incar=None):
    """Writes out the INCAR file.

    Boolean values are written as .TRUE./.FALSE.
    integers/floats and strings are written out as is
    lists/tuples are written out as space separated values/

    """

    if incar is None:
        incar = os.path.join(self.directory, 'INCAR')

    incar_keys = list(set(self.parameters) - set(self.special_kwargs))
    d = {key: self.parameters[key] for key in incar_keys}

    with open(incar, 'w') as f:
        f.write('INCAR created by Atomic Simulation Environment\n')
        for key, val in d.iteritems():
            key = ' ' + key.upper()
            if val is None:
                # Do not write out None values
                # It is how we delete tags
                pass
            elif isinstance(val, bool):
                s = '.TRUE.' if val else '.FALSE.'
                f.write('{} = {}\n'.format(key, s))
            elif isinstance(val, list) or isinstance(val, tuple):
                s = ' '.join([str(x) for x in val])
                f.write('{} = {}\n'.format(key, s))
            else:
                f.write('{} = {}\n'.format(key, val))


@monkeypatch_class(vasp.Vasp)
def write_kpoints(self, fname=None):
    """Write out the KPOINTS file.

    The KPOINTS file format is as follows:

    line 1: a comment
    line 2: number of kpoints
        n <= 0   Automatic kpoint generation
        n > 0    explicit number of kpoints
    line 3: kpt format
        if n > 0:
            C,c,K,k = cartesian coordinates
            anything else = reciprocal coordinates
        if n <= 0
            M,m,G,g for Monkhorst-Pack or Gamma grid
            anything else is a special case
    line 4: if n <= 0, the Monkhorst-Pack grid
        if n > 0, then a line per kpoint
    line 5: if n <=0 it is the gamma shift

    After the kpts may be tetrahedra, but we do now support that for
    now.

    """
    if fname is None:
        fname = os.path.join(self.directory, 'KPOINTS')

    p = self.parameters

    kpts = p.get('kpts', None)  # this is a list, or None

    if kpts is None:
        NKPTS = None
    elif len(np.array(kpts).shape) == 1:
        NKPTS = 0  # automatic
    else:
        NKPTS = len(p['kpts'])

    # figure out the mode
    if NKPTS == 0 and not p.get('gamma', None):
        MODE = 'm'  # automatic monkhorst-pack
    elif NKPTS == 0 and p.get('gamma', None):
        MODE = 'g'  # automatic gamma monkhorst pack
    # we did not trigger automatic kpoints
    elif p.get('kpts_nintersections', None) is not None:
        MODE = 'l'
    elif p.get('reciprocal', None) is True:
        MODE = 'r'
    else:
        MODE = 'c'

    with open(fname, 'w') as f:
        # line 1 - comment
        f.write('KPOINTS created by Atomic Simulation Environment\n')
        # line 2 - number of kpts
        if MODE in ['c', 'k', 'm', 'g', 'r']:
            f.write('{}\n'.format(NKPTS))
        elif MODE in ['l']:  # line mode, default intersections is 10
            f.write('{}\n'.format(p.get('kpts_nintersections')))

        # line 3
        if MODE in ['m', 'g']:
            if MODE == 'm':
                f.write('Monkhorst-Pack\n')  # line 3
            elif MODE == 'g':
                f.write('Gamma\n')
        elif MODE in ['c', 'k']:
            f.write('Cartesian\n')
        elif MODE in ['l']:
            f.write('Line-mode\n')
        else:
            f.write('Reciprocal\n')

        # line 4
        if MODE in ['m', 'g']:
            f.write('{0} {1} {2}\n'.format(*p.get('kpts', (1, 1, 1))))
        elif MODE in ['c', 'k', 'r']:
            for n in range(NKPTS):
                # I assume you know to provide the weights
                f.write('{0} {1} {2} {3}\n'.format(*p['kpts'][n]))
        elif MODE in ['l']:
            if p.get('reciprocal', None) is False:
                f.write('Cartesian\n')
            else:
                f.write('Reciprocal\n')
            for n in range(NKPTS):
                f.write('{0} {1} {2}\n'.format(*p['kpts'][n]))

        # line 5 - only if we are in automatic mode
        if MODE in ['m', 'g']:
            if p.get('gamma', None) and len(p.get('gamma')) == 3:
                f.write('{0} {1} {2}\n'.format(*p['gamma']))
            else:
                f.write('0.0 0.0 0.0\n')


@monkeypatch_class(vasp.Vasp)
def write_potcar(self, fname=None):
    """Writes the POTCAR file.

    POTCARs are expected in $VASP_PP_PATH.

    """
    if fname is None:
        fname = os.path.join(self.directory, 'POTCAR')

    with open(fname, 'wb') as potfile:
        for _, pfile, _ in self.ppp_list:
            pfile = os.path.join(os.environ['VASP_PP_PATH'], pfile)
            with open(pfile) as f:
                potfile.write(f.read())
