"""Reader functions for the Vasp calculator."""

import os
import re
import numpy as np
import vasp
from vasp import log
from ase.calculators.calculator import Parameters

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
def read_metadata(self, fname=None):
    """Read the METADATA file.

    Sets self.metadata to the results.

    """
    if fname is None:
        fname = os.path.join(self.directory, 'METADATA')

    self.metadata = {}
    if os.path.exists(fname):
        with open(fname) as f:
            import json
            metadata = json.loads(f.read())
            self.metadata = metadata
            log.debug('set metadata to: {}'.format(metadata))
    else:
        log.debug('{} did not exist.'.format(fname))


@monkeypatch_class(vasp.Vasp)
def read(self):
    """Read the files in a calculation.

    Returns the atoms and a Parameters dictionary.
    """

    params = Parameters()

    params.update(self.read_incar())
    params.update(self.read_potcar())
    params.update(self.read_kpoints())

    self.read_metadata()

    # We have to figure out the xc that was used based on the
    # Parameter keys.  We sort the possible xc dictionaries so the
    # ones with the largest number of keys are compared first. This is
    # to avoid false matches of xc's with smaller number of equal
    # keys.
    xc_keys = sorted(vasp.Vasp.xc_defaults,
                     key=lambda k: len(vasp.Vasp.xc_defaults[k]),
                     reverse=True)

    for ex in xc_keys:
        pd = {k: params.get(k, None) for k in vasp.Vasp.xc_defaults[ex]}
        if pd == vasp.Vasp.xc_defaults[ex]:
            params['xc'] = ex
            break

    import ase.io
    contcar = os.path.join(self.directory, 'CONTCAR')
    poscar = os.path.join(self.directory, 'POSCAR')

    if os.path.exists(contcar):
        atoms = ase.io.read(contcar)[self.metadata['resort']]
    else:
        atoms = ase.io.read(poscar)[self.metadata['resort']]
    return (atoms, params)
