"""Calculate bandstructure diagrams in jasp"""
import vasp
from monkeypatch import monkeypatch_class

import os
import numpy as np
import matplotlib.pyplot as plt
from ase.dft import DOS


@monkeypatch_class(vasp.Vasp)
def get_bandstructure(self,
                      kpts_path=None,
                      kpts_nintersections=10):
    """Calculate band structure along :param kpts_path:
    :param list kpts_path: list of tuples of (label, k-point) to
      calculate path on.
    :param int kpts_nintersections: is the number of points between
      points in band structures. More makes the bands smoother.

    returns (npoints, band_energies, fighandle)

    """
    self.stop_if(self.potential_energy is None)

    kpts = [k[1] for k in kpts_path]
    labels = [k[0] for k in kpts_path]

    dos = DOS(self, width=0.2)
    d = dos.get_dos()
    e = dos.get_energies()

    ef = self.get_fermi_level()

    # run in non-selfconsistent directory

    wd = os.path.join(self.directory, 'bandstructure')
    self.clone(wd)

    calc = vasp.Vasp(wd)
    calc.set(kpts=kpts,
             kpts_nintersections=kpts_nintersections,
             reciprocal=True,
             nsw=0,  # no ionic updates required
             isif=None,
             ibrion=None,
             icharg=11)
    calc.update()

    calc.stop_if(calc.get_potential_energy() is None)

    fig = plt.figure()
    with open(os.path.join(wd, 'EIGENVAL')) as f:
        # skip 5 lines
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        unknown, npoints, nbands = [int(x) for x in f.readline().split()]

        f.readline()  # skip line

        band_energies = [[] for i in range(nbands)]

        for i in range(npoints):
            x, y, z, weight = [float(x) for x in f.readline().split()]

            for j in range(nbands):
                fields = f.readline().split()
                id, energy = int(fields[0]), float(fields[1])
                band_energies[id - 1].append(energy)
            f.readline()  # skip line

    ax1 = plt.subplot(121)
    for i in range(nbands):
        plt.plot(range(npoints), np.array(band_energies[i]) - ef)

    ax = plt.gca()
    ax.set_xticks([])  # no tick marks
    plt.xlabel('k-vector')
    plt.ylabel('Energy (eV)')

    nticks = len(labels) / 2 + 1
    ax.set_xticks(np.linspace(0, npoints, nticks))
    L = []
    L.append(labels[0])
    for i in range(2, len(labels)):
        if i % 2 == 0:
            L.append(labels[i])
        else:
            pass
    L.append(labels[-1])
    ax.set_xticklabels(L)
    plt.axhline(0, c='r')

    plt.subplot(122, sharey=ax1)
    plt.plot(d, e)
    plt.axhline(0, c='r')
    plt.ylabel('energy (eV)')
    plt.xlabel('DOS')

    plt.subplots_adjust(wspace=0.26)
    plt.show()
    return (npoints, band_energies, fig)
