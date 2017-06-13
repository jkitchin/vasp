'''module for vibrational calculations in jasp'''

import os
import numpy as np
import vasp
from .vasp import log, Vasp
from .monkeypatch import monkeypatch_class

from ase.data import atomic_masses
from ase.io import read


@monkeypatch_class(Vasp)
def get_vibrational_modes(self,
                          mode=None,
                          massweighted=False,
                          show=False,
                          npoints=30,
                          amplitude=0.5):

    """Read the OUTCAR and get the eigenvectors. Return value depends
    on the arguments.

    mode= None returns all modes
    mode= 2 returns mode 2
    mode=[1, 2] returns modes 1 and 2

    massweighted = True returns sqrt(mass) weighted
    eigenvectors. E.g. M * evectors * M

    show=True makes a trajectory that can be visualized
    npoints = number of points in the trajectory
    amplitude = magnitude of the vibrations

    some special cases to handle:
    ibrion=5 + selective dynamics
       may lead to unexpected number of modes

    if nwrite=3, there will be a sqrt(mass) weighted vectors
    and two sets of vectors.

    I am not sure if these eigenvectors are mass-weighted. And I am
    not sure if the order of the eigenvectors in OUTCAR is the same as
    the atoms.

    Note: it seems like it might be much easier to get this out of
    vasprun.xml
    """
    self.update()

    atoms = self.get_atoms()

    if hasattr(atoms, 'constraints') and self.parameters['ibrion'] == 5:
        # count how many modes to get.
        NMODES = 0
        f = open(os.path.join(self.directory, 'OUTCAR'))
        for line in f:
            if ('f' in line and 'THz' in line and 'cm-1' in line):
                NMODES += 1
        f.close()
    else:
        NMODES = 3 * len(atoms)

    frequencies, eigenvectors = [], []

    # now we find where the data starts. I think the unweighted
    # vectors always come first. if nwrite=3, then there are
    # sqrt(mass) weighted vectors that follow this section

    f = open(os.path.join(self.directory, 'OUTCAR'), 'r')
    while True:
        line = f.readline()
        if line.startswith(' Eigenvectors and eigenvalues'
                           ' of the dynamical matrix'):
            break
    f.readline()   # skip ------
    f.readline()   # skip two blank lines
    f.readline()

    for i in range(NMODES):
        freqline = f.readline()
        fields = freqline.split()

        if 'f/i=' in freqline:  # imaginary frequency
            frequencies.append(complex(float(fields[-2]) * 0.001, 0j))
        else:
            frequencies.append(float(fields[-2]) * 0.001)
        #        X         Y         Z           dx          dy          dz
        f.readline()
        thismode = []
        for i in range(len(atoms)):
            line = f.readline().strip()
            X, Y, Z, dx, dy, dz = [float(x) for x in line.split()]
            thismode.append(np.array([dx, dy, dz]))
        f.readline()  # blank line

        thismode = np.array(thismode)
        # now we need to resort the vectors in this mode so they match
        # the atoms order
        thismode = thismode[self.resort]

        if massweighted:
            # construct M
            numbers = [a.get('number') for a in atoms]
            M = []
            for i in range(len(atoms)):
                for j in range(3):
                    an = numbers[i]
                    M.append(1. / np.sqrt(atomic_masses[an]))
            M = np.array(M)
            M = np.diag(M)  # diagonal array

            thismode = np.dot(M, thismode.flat)

            thismode = thismode.reshape((len(atoms), 3))
        # renormalize the mode
        mag = np.linalg.norm(thismode)
        thismode /= mag

        eigenvectors.append(thismode)
    f.close()

    eigenvectors = np.array(eigenvectors)

    if mode is None:
        retval = (frequencies, eigenvectors)
    else:
        retval = (frequencies[mode], eigenvectors[mode])

    if show:
        from ase.visualize import view
        if mode is None:
            mode = [0]
        elif not isinstance(mode, list):
            mode = [mode]  # make a list for next code

        # symmetric path from -1 to 1 to -1
        X = np.append(np.linspace(0, 1, npoints / 3),
                      np.linspace(1, -1, npoints / 3))
        X = np.append(X,
                      np.linspace(-1, 0, npoints / 3))
        X *= amplitude

        for m in mode:
            traj = []
            for i, x in enumerate(X):
                a = atoms.copy()
                a.positions += x * eigenvectors[m]
                traj += [a]

            view(traj)
    return retval


@monkeypatch_class(Vasp)
def get_vibrational_frequencies(self):
    """Returns an array of frequencies in wavenumbers.

    You should have run the calculation already. This function does not
    run a calculation.
    """
    self.update()
    atoms = self.get_atoms()
    N = len(atoms)

    frequencies = []

    f = open(os.path.join(self.directory, 'OUTCAR'), 'r')
    while True:
        line = f.readline()
        if line.startswith(' Eigenvectors and eigenvalues'
                           ' of the dynamical matrix'):
            break
    f.readline()  # skip ------
    f.readline()  # skip two blank lines
    f.readline()
    for i in range(3 * N):
        # the next line contains the frequencies
        line = f.readline()
        fields = line.split()

        if 'f/i=' in line:  # imaginary frequency
            # frequency in wave-numbers
            frequencies.append(complex(float(fields[6]), 0j))
        else:
            frequencies.append(float(fields[7]))
        # now skip 1 one line, a line for each atom, and a blank line
        for j in range(1 + N + 1):
            f.readline()  # skip the next few lines
    f.close()
    return frequencies


@monkeypatch_class(Vasp)
def get_infrared_intensities(self):
    """Calculate infrared intensities of vibrational modes.

    Returns an array of normalized intensities for each vibrational
    mode. You should have run the vibrational calculation already. This
    function does not run it for you.

    python translation of # A utility for calculating the vibrational
    intensities from VASP output (OUTCAR) # (C) David Karhanek,
    2011-03-25, ICIQ Tarragona, Spain (www.iciq.es)
    http://homepage.univie.ac.at/david.karhanek/downloads.html#Entry02
    """
    assert self.parameters.get('lepsilon', None) is True
    assert self.parameters.get('nwrite', 0) == 3
    assert self.parameters.get('ibrion', 0) == 7

    self.update()

    atoms = read(os.path.join(self.directory, 'POSCAR'), format='vasp')
    NIONS = len(atoms)
    BORN_NROWS = NIONS * 4 + 1

    with open(os.path.join(self.directory, 'OUTCAR'), 'r') as f:
        alltext = f.read()
        f.seek(0)
        alllines = f.readlines()
        f.close()

    if 'BORN' not in alltext:
        raise Exception('Born effective charges missing. '
                        'Did you use IBRION=7 or 8?')

    if 'Eigenvectors after division by SQRT(mass)' not in alltext:
        raise Exception('You must rerun with NWRITE=3 to get '
                        'sqrt(mass) weighted eigenvectors')

    # get the Born charges
    for i, line in enumerate(alllines):
        if 'BORN EFFECTIVE CHARGES' in line:
            break

    BORN_MATRICES = []
    i += 2  # skip a line
    for j in range(NIONS):
        BM = []
        i += 1  # skips the ion count line
        for k in range(3):
            line = alllines[i]
            fields = line.split()
            BM.append([float(x) for x in fields[1:4]])
            i += 1  # advance a line
        BORN_MATRICES.append(BM)

    BORN_MATRICES = np.array(BORN_MATRICES)

    # Get the eigenvectors and eigenvalues.  maybe I can replace this
    # code with my other code. for now I just reproduce the count
    # number of vibs. this gets the number from outcar. it seems like
    # it should be known in advance unless constraints make it hard to
    # tell.

    # the next code in the shell script just copies code to eigenvectors.txt
    for i, line in enumerate(alllines):
        if 'Eigenvectors after division by SQRT(mass)' in line:
            break

    EIG_NVIBS = 0
    for line in alllines[i:]:
        if ('f' in line and
            'THz' in line and
            'cm-1' in line):
            EIG_NVIBS += 1

    EIG_NIONS = BORN_NROWS
    # I guess this counts blank rows and non-data rows
    # EIG_NROWS = (EIG_NIONS + 3) * EIG_NVIBS + 3

    # i is where the data starts
    i += 6

    EIGENVALUES = []
    EIGENVECTORS = []
    for j in range(EIG_NVIBS):
        mode = []
        EIGENVALUES.append(alllines[i])  # frequencies are here

        i += 1  # skip the frequency line
        i += 1  # skip the xyz line
        for k in range(3):
            fields = [float(x) for x in alllines[i].split()]
            mode.append(fields[3:])
            i += 1
        EIGENVECTORS.append(mode)
        i += 1  # skip blank line

    EIGENVECTORS = np.array(EIGENVECTORS)

    # now we are ready to compute intensities. see
    # http://othes.univie.ac.at/10117/1/2010-05-05_0547640.pdf, page
    # 21.

    # I(\omega) = \sum_{\alpha=1}^3 |
    # \sum_{l=1}^M \sum_{\beta=1}^3 Z_{\alpha\beta}(l)e_{\beta}(l)|^2

    # omega is the vibrational mode
    # alpha, beta are the cartesian polarizations
    # l is the atom number
    # e_beta is the eigenvector of the mode

    intensities = []

    for mode in range(len(EIGENVECTORS)):
        S = 0  # This is the triple sum
        for alpha in [0, 1, 2]:
            s = 0
            for l in [0, 1, 2]:  # this is the atom number
                for beta in [0, 1, 2]:
                    e = EIGENVECTORS[mode][l]
                    Zab = BORN_MATRICES[l][alpha][beta]

                    s += Zab * e[beta]
            S += s ** 2
        intensities.append(S)

    intensities = np.array(intensities) / max(intensities)
    return intensities
