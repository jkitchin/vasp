import os
from subprocess import Popen, PIPE
import vasp
from monkeypatch import monkeypatch_class


@monkeypatch_class(vasp.Vasp)
def chgsum(self):
    """Uses the chgsum.pl utility to sum over the AECCAR0 and AECCAR2 files."""
    cmdlist = ['chgsum.pl', 'AECCAR0', 'AECCAR2']
    p = Popen(cmdlist, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    if out == '' or err != '':
        raise Exception('Cannot perform chgsum:\n\n{0}'.format(err))


@monkeypatch_class(vasp.Vasp)
def bader(self, cmd=None, ref=False, verbose=False, overwrite=False):
    """Performs bader analysis for a calculation.
    Follows defaults unless full shell command is specified
    Does not overwrite existing files if overwrite=False
    If ref = True, tries to reference the charge density to
    the sum of AECCAR0 and AECCAR2
    Requires the bader.pl (and chgsum.pl) script to be in the system PATH
    """

    if 'ACF.dat' in os.listdir('./') and not overwrite:
        return

    if cmd is None:
        if ref:
            self.chgsum()
            cmdlist = ['bader', 'CHGCAR', '-ref', 'CHGCAR_sum']
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
        print('Bader completed for {0}'.format(self.vaspdir))
