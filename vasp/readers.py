"""Reader functions for the Vasp calculator."""

import os
import re
import numpy as np
import vasp
from vasp import log
from ase.calculators.calculator import Parameters
import exceptions
from monkeypatch import monkeypatch_class


def isfloat(s):
    """Return if s is a float.

    We check if it is not an integer first, then try to make it a float.

    """
    if re.match('^[-+]?\d+$', s):
        return False
    try:
        float(s)
        return True
    except:
        return False


@monkeypatch_class(vasp.Vasp)
def read_incar(self, fname=None):
    """Read fname (defaults to INCAR).

    Returns a Parameters dictionary from the INCAR.

    This only reads simple INCAR files, e.g. one tag per line, and
    with no comments in the line. There is no knowledge of any Vasp
    keywords in this, and the values are converted to Python types by
    some simple rules.

    """

    if fname is None:
        fname = self.incar

    params = Parameters()

    with open(fname) as f:
        lines = f.readlines()

    # The first line is a comment
    for line in lines[1:]:
        line = line.strip()
        if ";" in line:
            raise Exception('; found. that is not supported.')
        if '#' in line:
            raise Exception('# found. that is not supported.')
        if line == '':
            continue

        key, val = line.split('=')
        key = key.strip().lower()
        val = val.strip()
        # now we need some logic.
        if val == '.TRUE.':
            val = True
        elif val == '.FALSE.':
            val = False
        # Match integers by a regexp that includes signs
        # val.isdigit() does not get negative integers right.
        elif re.match('^[-+]?\d+$', val):
            val = int(val)
        elif isfloat(val):
            val = float(val)
        elif len(val.split(' ')) > 1:
            # this is some kind of list separated by spaces
            val = val.split(' ')
            val = [int(x) if re.match('^[-+]?\d+$', x)
                   else float(x) for x in val]
        else:
            # I guess we have a string here.
            pass

        params[key] = val

    return params


@monkeypatch_class(vasp.Vasp)
def read_kpoints(self, fname=None):
    """Read KPOINTS file.

    Returns a Parameters object of kpoint tags.

    """

    if fname is None:
        fname = self.kpoints

    with open(fname) as f:
        lines = f.readlines()

    params = Parameters()

    # first line is a comment
    # second line is the number of kpoints or 0 for automatic kpoints
    nkpts = int(lines[1].strip())

    # third line you have to specify whether the coordinates are given
    # in cartesian or reciprocal coordinates if nkpts is greater than
    # zero. Only the first character of the third line is
    # significant. The only key characters recognized by VASP are 'C',
    # 'c', 'K' or 'k' for switching to cartesian coordinates, any
    # other character will switch to reciprocal coordinates.
    #
    # if nkpts = 0 then the third line will start with m or g for
    # Monkhorst-Pack and Gamma. if it does not start with m or g, an
    # alternative mode is entered that we do not support yet.

    ktype = lines[2].split()[0].lower()[0]
    if nkpts <= 0:
        # automatic mode
        if ktype not in ['g', 'm']:
            raise NotImplementedError('Only Monkhorst-Pack and '
                                      'gamma centered grid supported '
                                      'for restart.')
        if ktype == 'g':
            line5 = np.array([float(lines[4].split()[i]) for i in range(3)])
            if (line5 == np.array([0.0, 0.0, 0.0])).all():
                params['gamma'] = True
            else:
                params['gamma'] = line5

        kpts = [int(lines[3].split()[i]) for i in range(3)]
        params['kpts'] = kpts
    elif nkpts > 0:
        # list of kpts provided. Technically c,k are supported and
        # anything else means reciprocal coordinates.
        if ktype in ['c', 'k', 'r']:
            kpts = []
            for i in range(3, 3 + nkpts):
                # the kpts also have a weight attached to them
                kpts.append([float(lines[i].split()[j])
                             for j in range(4)])
            params['kpts'] = kpts
        # you may also be in line-mode
        elif ktype in ['l']:
            if lines[3][0].lower() == 'r':
                params['reciprocal'] = True

            params['kpts_nintersections'] = nkpts

            kpts = []
            for i in range(4, len(lines)):
                if lines[i] == '':
                    continue
                else:
                    kpts.append([float(lines[i].split()[j])
                                 for j in range(3)])
        else:
            raise NotImplementedError('ktype = %s' % lines[2])

    if ktype == 'r':
        params['reciprocal'] = True

    params['kpts'] = kpts
    return params


@monkeypatch_class(vasp.Vasp)
def read_potcar(self, fname=None):
    """Read the POTCAR file to get the pp and setups.

    Returns a Parameters dictionary of pp and setups.

    """

    if fname is None:
        fname = self.potcar

    params = Parameters()

    potcars = []
    with open(fname) as f:
        lines = f.readlines()

    # first potcar
    potcars += [lines[0].strip()]

    for i, line in enumerate(lines):
        if 'LEXCH  = PE' in line:
            params['pp'] = 'PBE'
        elif 'LEXCH  = CA' in line:
            params['pp'] = 'LDA'
        elif 'LEXCH  = 91' in line:
            params['pp'] = 'GGA'

        if 'End of Dataset' in line and i != len(lines) - 1:
            potcars += [lines[i + 1].strip()]

    potcars = [(x[0], x[1], x[2]) for x in
               [potcar.split() for potcar in potcars]]

    special_setups = []
    for xc, sym, date in potcars:
        if '_' in sym:  # we have a special setup
            symbol, setup = sym.split('_')
            special_setups += [[symbol,  '_' + setup]]

    if special_setups:
        params['setups'] = special_setups

    return params


@monkeypatch_class(vasp.Vasp)
def read(self, restart=None):
    """Read the files in a calculation if they exist.

    restart is ignored, but part of the signature for ase. I am not
    sure what we could use it for.

    sets self.parameters and atoms.

    """

    self.neb = None
    # NEB is special and handled separately
    if self.get_state() == vasp.Vasp.NEB:
        self.read_neb()
        return

    # Else read a regular calculation. we start with reading stuff
    # that is independent of the calculation state.
    self.parameters = Parameters()

    if os.path.exists(self.incar):
        self.parameters.update(self.read_incar())
    if os.path.exists(self.potcar):
        self.parameters.update(self.read_potcar())
    if os.path.exists(self.kpoints):
        self.parameters.update(self.read_kpoints())

    # We have to figure out the xc that was used based on the
    # Parameter keys.  We sort the possible xc dictionaries so the
    # ones with the largest number of keys are compared first. This is
    # to avoid false matches of xc's with smaller number of equal
    # keys.
    xc_keys = sorted(vasp.Vasp.xc_defaults,
                     key=lambda k: len(vasp.Vasp.xc_defaults[k]),
                     reverse=True)

    for ex in xc_keys:
        pd = {k: self.parameters.get(k, None)
              for k in vasp.Vasp.xc_defaults[ex]}
        if pd == vasp.Vasp.xc_defaults[ex]:
            self.parameters['xc'] = ex
            break

    # reconstruct ldau_luj. special setups might break this.
    if 'ldauu' in self.parameters:
        ldaul = self.parameters['ldaul']
        ldauj = self.parameters['ldauj']
        ldauu = self.parameters['ldauu']

        with open(self.potcar) as f:
            lines = f.readlines()

        # symbols are in the first line of each potcar
        symbols = [lines[0].split()[1]]
        for i, line in enumerate(lines):
            if 'End of Dataset' in line and i != len(lines) - 1:
                symbols += [lines[i + 1].split()[1]]

        ldau_luj = {}
        for sym, l, j, u in zip(symbols, ldaul, ldauj, ldauu):
            ldau_luj[sym] = {'L': l, 'U': u, 'J': j}

        self.parameters['ldau_luj'] = ldau_luj

    # Now for the atoms. This does depend on the state. self.resort
    # needs to be a list for shuffling constraints if they exist.
    self.resort = self.get_db('resort')
    if self.resort is not None:
        self.resort = list(self.resort)

    import ase.io
    contcar = os.path.join(self.directory, 'CONTCAR')
    empty_contcar = False
    if os.path.exists(contcar):
        # make sure the contcar is not empty
        with open(contcar) as f:
            if f.read() == '':
                empty_contcar = True

    poscar = os.path.join(self.directory, 'POSCAR')

    if os.path.exists(contcar) and not empty_contcar:
        atoms = ase.io.read(contcar)
    elif os.path.exists(poscar):
        atoms = ase.io.read(poscar)
    else:
        atoms = None

    if atoms is not None:
        atoms = atoms[self.resort]
        self.sort_atoms(atoms)

    self.read_results()


@monkeypatch_class(vasp.Vasp)
def read_results(self):
    """Read energy, forces, stress, magmom and magmoms from output file.

    Other quantities will be read by other functions. This depends on
    state.

    """
    state = self.get_state()
    if state == vasp.Vasp.NEB:
        # This is handled in self.read()
        return

    if state != vasp.Vasp.FINISHED:
        self.results = {}
    else:
        # regular calculation that is finished
        from ase.io.vasp import read_vasp_xml
        if not os.path.exists(os.path.join(self.directory,
                                           'vasprun.xml')):
            exc = 'No vasprun.xml in {}'.format(self.directory)
            raise exceptions.VaspNotFinished(exc)

        atoms = read_vasp_xml(os.path.join(self.directory,
                                           'vasprun.xml')).next()

        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()  # needs to be resorted
        stress = atoms.get_stress()

        resort = self.get_db('resort')
        if self.atoms is None:
            atoms = atoms[resort]
            self.sort_atoms(atoms)
            self.atoms.set_calculator(self)
        else:
            # update the atoms
            self.atoms.positions = atoms.positions[resort]
            self.atoms.cell = atoms.cell

        self.results['energy'] = energy
        self.results['forces'] = forces[self.resort]
        self.results['stress'] = stress
        self.results['dipole'] = None
        self.results['charges'] = np.array([None for atom in self.atoms])

        magnetic_moment = 0
        magnetic_moments = np.zeros(len(atoms))
        if self.parameters.get('ispin', 0) == 2:
            lines = open(os.path.join(self.directory, 'OUTCAR'),
                         'r').readlines()
            for n, line in enumerate(lines):
                if line.rfind('number of electron  ') > -1:
                    magnetic_moment = float(line.split()[-1])

                if line.rfind('magnetization (x)') > -1:
                    for m in range(len(atoms)):
                        magnetic_moments[m] = float(lines[n + m + 4].split()[4])

        self.results['magmom'] = magnetic_moment
        self.results['magmoms'] = np.array(magnetic_moments)[self.resort]


@monkeypatch_class(vasp.Vasp)
def read_neb(self):
    """Read an NEB calculator."""
    import ase
    import glob
    atoms = []
    atoms += [ase.io.read('{}/00/POSCAR'.format(self.directory))]
    for p in glob.glob('{}/0[0-9]/CONTCAR'.format(self.directory)):
        atoms += [ase.io.read(p)]
    atoms += [ase.io.read('{}/0{}/POSCAR'.format(self.directory,
                                                 len(atoms)))]
    self.neb = atoms
    self.parameters = {}
    self.set(images=(len(atoms) - 2))
    self.atoms = atoms[0].copy()

    if os.path.exists(self.incar):
        self.parameters.update(self.read_incar())
    if os.path.exists(self.potcar):
        self.parameters.update(self.read_potcar())
    if os.path.exists(self.kpoints):
        self.parameters.update(self.read_kpoints())


    # Update the xc functional
    xc_keys = sorted(vasp.Vasp.xc_defaults,
                     key=lambda k: len(vasp.Vasp.xc_defaults[k]),
                     reverse=True)

    for ex in xc_keys:
        pd = {k: self.parameters.get(k, None)
              for k in vasp.Vasp.xc_defaults[ex]}
        if pd == vasp.Vasp.xc_defaults[ex]:
            self.parameters['xc'] = ex
            break
