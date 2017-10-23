"""
code for running NEB calculations in Vasp

here is typical code to set up the band:

calc = Vasp('../surfaces/Pt-slab-O-fcc')
initial_atoms = calc.get_atoms()

calc = Vasp('../surfaces/Pt-slab-O-hcp')
final_atoms = calc.get_atoms()

images = [initial_atoms]
images += [initial_atoms.copy() for i in range(3)]
images += [final_atoms]

neb = NEB(images)
# Interpolate linearly the positions of the three middle images:
neb.interpolate()

calc = Vasp('O-diffusion',
            ibrion=2,
            nsw=50,
            images=5,  # initial + nimages + final
            spring=-5,
            atoms=images)
images, energies = calc.get_neb()

The spring tag triggers the setup of an NEB calculation for Vasp.

"""

import os
import numpy as np

from ase.io import read
from ase.io.vasp import write_vasp
from ase.calculators.vasp import Vasp

import vasp
import vasp.exceptions
from .vasp import log, Vasp
from .monkeypatch import monkeypatch_class
from .vasprc import VASPRC


@monkeypatch_class(Vasp)
def get_neb(self, npi=1):
    """Returns images, energies if available or runs the job.

    npi = cores per image for running the calculations. Default=1

    show: if True show an NEB plot
    """
    if self.in_queue():
        return self.neb, [None for a in self.neb]

    calc_required = False

    # check for OUTCAR in each image dir
    for i in range(1, len(self.neb) - 1):
        wf = '{0}/OUTCAR'.format(str(i).zfill(2))
        wf = os.path.join(self.directory, wf)
        if not os.path.exists(wf):
            calc_required = True
            break
        else:
            # there was an OUTCAR, now we need to check for
            # convergence.
            done = False
            with open(wf) as f:
                for line in f:
                    if ('reached required accuracy - stopping structural'
                        ' energy minimisation') in line:
                        done = True
                        break
            if not done:
                calc_required = True
                break

    if calc_required:
        log.debug('NEB calculation required')
        # this creates the directories and files if needed.  write out
        # all the images, including initial and final
        if not os.path.isdir(self.directory):
            os.makedirs(self.directory)

        self.set(images=len(self.neb) - 2)

        self.write_incar()
        log.debug('Wrote incar')
        self.write_kpoints()
        log.debug('Wrote kpoints')
        self.write_potcar()
        log.debug('Wrote potcar')

        # This is hanging
        self.write_db();
        log.debug('Wrote db')

        for i, atoms in enumerate(self.neb):
            # zero-padded directory name
            image_dir = os.path.join(self.directory, str(i).zfill(2))
            if not os.path.isdir(image_dir):
                # create if needed.
                log.debug('Creating {}'.format(image_dir))
                os.makedirs(image_dir)
                write_vasp('{0}/POSCAR'.format(image_dir),
                           atoms[self.resort],
                           symbol_count=self.symbol_count)

        # The first and last images need to have real calculators on
        # them so we can write out a DB entry. We need this so we can
        # get the energies on the end-points. Otherwise, there doesn't
        # seem to be a way to do that short of cloning the whole
        # calculation into the end-point directories.

        log.debug('Writing initial state db')
        self.write_db(os.path.join(self.directory,
                                   '00/DB.db'),
                      self.neb[0])

        log.debug('Writing final state db')
        self.write_db(os.path.join(self.directory,
                                   '{:02}/DB.db'.format(len(self.neb) - 1)),
                      self.neb[-1])

        VASPRC['queue.ppn'] = npi * (len(self.neb) - 2)
        log.debug('Running on %i cores', VASPRC['queue.ppn'])

        self.calculate()  # this will raise VaspSubmitted
        return self.neb, [None for a in self.neb]

    #############################################
    # now we are just retrieving results
    energies = []
    import ase.io
    atoms0 = ase.io.read(os.path.join(self.directory,
                                      '00',
                                      'DB.db'))
    energies += [atoms0.get_potential_energy()]

    for i in range(1, len(self.neb) - 1):
        atoms = ase.io.read(os.path.join(self.directory,
                                         str(i).zfill(2),
                                         'CONTCAR'))[self.resort]
        self.neb[i].positions = atoms.positions
        self.neb[i].cell = atoms.cell

        energy = None
        with open(os.path.join(self.directory,
                               str(i).zfill(2),
                               'OUTCAR')) as f:
            for line in f:
                if 'free energy    TOTEN  =' in line:
                    energy = float(line.split()[4])

        energies += [energy]

    fname = os.path.join(self.directory,
                         '0{}/DB.db'.format(len(self.neb) - 1))
    atoms_end = ase.io.read(fname)
    energies += [atoms_end.get_potential_energy()]

    energies = np.array(energies)
    energies -= energies[0]

    return (self.neb, np.array(energies))


@monkeypatch_class(Vasp)
def plot_neb(self, show=True):
    """Return a list of the energies and atoms objects for each image in

    the band.

    by default shows the plot figure
    """
    images, energies = self.get_neb()
    # add fitted line to band energies. we make a cubic spline
    # interpolating function of the negative energy so we can find the
    # minimum which corresponds to the barrier
    from scipy.interpolate import interp1d
    from scipy.optimize import fmin
    f = interp1d(list(range(len(energies))),
                 -energies,
                 kind='cubic', bounds_error=False)
    x0 = len(energies) / 2.  # guess barrier is at half way
    xmax = fmin(f, x0)

    xfit = np.linspace(0, len(energies) - 1)
    bandfit = -f(xfit)

    import matplotlib.pyplot as plt
    p = plt.plot(energies - energies[0], 'bo ', label='images')
    plt.plot(xfit, bandfit, 'r-', label='fit')
    plt.plot(xmax, -f(xmax), '* ', label='max')
    plt.xlabel('Image')
    plt.ylabel('Energy (eV)')
    s = ['$\Delta E$ = {0:1.3f} eV'.format(float(energies[-1] -
                                                 energies[0])),
         '$E^\ddag$ = {0:1.3f} eV'.format(float(-f(xmax)))]

    plt.title('\n'.join(s))
    plt.legend(loc='best', numpoints=1)
    if show:
        from ase.calculators.singlepoint import SinglePointCalculator
        from ase.visualize import view
        # It seems there might be some info on the atoms that causes
        # an error here. Making a copy seems to get rid of the
        # issue. Hacky.
        tatoms = [x.copy() for x in images]
        for i, x in enumerate(tatoms):
            x.set_calculator(SinglePointCalculator(x, energy=energies[i]))
        view(tatoms)
        plt.show()
    return p
