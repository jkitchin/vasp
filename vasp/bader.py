import os
import numpy as np
from subprocess import Popen, PIPE
import vasp
from POTCAR import get_ZVAL
from monkeypatch import monkeypatch_class


@monkeypatch_class(vasp.Vasp)
def chgsum(self):
    """Uses the chgsum.pl utility to sum over the AECCAR0 and AECCAR2 files."""
    cwd = os.getcwd()
    try:
        os.chdir(self.directory)
        cmdlist = ['chgsum.pl',
                   'AECCAR0',
                   'AECCAR2']
        p = Popen(cmdlist, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        if out == '' or err != '':
            raise Exception('Cannot perform chgsum:\n\n{0}'.format(err))
    finally:
        os.chdir(cwd)


@monkeypatch_class(vasp.Vasp)
def bader(self, cmd=None, ref=False, verbose=False, overwrite=False):
    """Performs bader analysis for a calculation.
    Follows defaults unless full shell command is specified
    Does not overwrite existing files if overwrite=False
    If ref = True, tries to reference the charge density to
    the sum of AECCAR0 and AECCAR2
    Requires the bader.pl (and chgsum.pl) script to be in the system PATH
    """
    cwd = os.getcwd()
    try:
        os.chdir(self.directory)

        if 'ACF.dat' in os.listdir(".") and not overwrite:
            self._get_calculated_charges()
            return

        if cmd is None:
            if ref:
                self.chgsum()
                cmdlist = ['bader',
                           'CHGCAR',
                           '-ref',
                           'CHGCAR_sum']
            else:
                cmdlist = ['bader', 'CHGCAR']
        elif type(cmd) is str:
            cmdlist = cmd.split()
        elif type(cmd) is list:
            cmdlist = cmd

        p = Popen(cmdlist, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        if out == '' or err != '':
            raise Exception('Cannot perform Bader:\n\n{0}'.format(err))
        elif verbose:
            print('Bader completed for {0}'.format(self.directory))

        # Now store the calculated charges
        self._get_calculated_charges()

    finally:
        os.chdir(cwd)


@monkeypatch_class(vasp.Vasp)
def _get_calculated_charges(self,
                            atoms=None,
                            fileobj='ACF.dat',
                            displacement=1e-4):

    """Calculate the charges from the fileobj.
    This is a modified version of the attach_charges function in
    ase.io.bader to work better with VASP.
    Does not require the atom positions to be in Bohr and references
    the charge to the ZVAL in the POTCAR
    """

    if fileobj is None:
        fileobj = os.path.join(self.directory, 'ACF.dat')

    if isinstance(fileobj, str):
            fileobj = open(fileobj)
            f_open = True

    if atoms is None:
        atoms = self.get_atoms()

    # Get the sorting and resorting lists
    sort = self.sort
    resort = self.resort

    # First get a dictionary of ZVALS from the pseudopotentials
    LOP = self.get_pseudopotentials()
    ppp = os.environ['VASP_PP_PATH']

    zval = {}
    for sym, ppath, hash in LOP:
        fullpath = ppp + ppath
        z = get_ZVAL(fullpath)
        zval[sym] = z

    # Get sorted symbols and positions according to POSCAR and ACF.dat
    symbols = np.array(atoms.get_chemical_symbols())[sort]
    positions = atoms.get_positions()[sort]

    charges = []
    sep = '---------------'
    i = 0  # Counter for the lines
    k = 0  # Counter of sep
    assume6columns = False
    for line in fileobj:
        if line[0] == '\n':  # check if there is an empty line in the
            i -= 1           # head of ACF.dat file
        if i == 0:
            headings = line
            if 'BADER' in headings.split():
                j = headings.split().index('BADER')
            elif 'CHARGE' in headings.split():
                j = headings.split().index('CHARGE')
            else:
                print('Can\'t find keyword "BADER" or "CHARGE".' \
                      + ' Assuming the ACF.dat file has 6 columns.')
                j = 4
                assume6columns = True
        if sep in line:  # Stop at last seperator line
            if k == 1:
                break
            k += 1
        if not i > 1:
            pass
        else:
            words = line.split()
            if assume6columns is True:
                if len(words) != 6:
                    raise IOError('Number of columns in ACF file incorrect!\n'
                                  'Check that Bader program version >= 0.25')

            sym = symbols[int(words[0]) - 1]
            charges.append(zval[sym] - float(words[j]))

            if displacement is not None:
                # check if the atom positions match
                xyz = np.array([float(w) for w in words[1:4]])
                assert np.linalg.norm(positions[int(words[0]) - 1] - xyz) < displacement
        i += 1

    if f_open:
        fileobj.close()

    # Now attach the resorted charges to the atom
    charges = np.array(charges)[resort]
    self._calculated_charges = charges
