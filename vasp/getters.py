import os
import numpy as np
from xml.etree import ElementTree
import vasp
from vasp import log
from monkeypatch import monkeypatch_class


@monkeypatch_class(vasp.Vasp)
def get_k_point_weights(self):
    """Return the k-point weights."""
    self.read_results()
    with open(os.path.join(self.directory,
                           'vasprun.xml')) as f:
        tree = ElementTree.parse(f)
        # each weight is in a <v>w</v> element in this varray
        return np.array([float(x.text) for x in
                         tree.find("kpoints/varray[@name='weights']")])


@monkeypatch_class(vasp.Vasp)
def get_number_of_spins(self):
    """Returns number of spins.
    1 if not spin-polarized
    2 if spin-polarized

    """
    if 'ispin' in self.parameters:
        return 2
    else:
        return 1


@monkeypatch_class(vasp.Vasp)
def get_eigenvalues(self, kpt=0, spin=1):
    """Return array of eigenvalues for kpt and spin."""
    log.debug('kpt={} spin={}'.format(kpt, spin))
    self.read_results()
    with open(os.path.join(self.directory,
                           'vasprun.xml')) as f:
        tree = ElementTree.parse(f)
        path = '/'.join(['calculation',
                         'eigenvalues',
                         'array',
                         'set',
                         "set[@comment='spin {}']",
                         "set[@comment='kpoint {}']"])
        path = path.format(spin + 1, kpt + 1)
        log.debug('path={}'.format(path))
        # Vasp seems to start at 1 not 0
        fields = tree.find(path)

        return np.array([float(x.text.split()[0]) for x in fields])


@monkeypatch_class(vasp.Vasp)
def get_fermi_level(self):
    """Return the Fermi level."""
    self.read_results()
    with open(os.path.join(self.directory,
                           'vasprun.xml')) as f:
        tree = ElementTree.parse(f)
        path = '/'.join(['calculation',
                         'dos',
                         "i[@name='efermi']"
                         ])
        return float(tree.find(path).text)


#@monkeypatch_class(vasp.Vasp)
#def get_magnetic_moment(self

@monkeypatch_class(vasp.Vasp)
def get_ados(self, atom_index, orbital, spin=1, efermi=None):
    """Return Atom projected DOS for atom index, orbital and spin.

    orbital: string ['s', 'p', 'd']

    If efermi is not None, use this value as 0.0.

    """
    self.read_results()

    with open(os.path.join(self.directory,
                           'vasprun.xml')) as f:
        tree = ElementTree.parse(f)

    path = "/".join(['calculation', 'dos',
                     'partial',
                     'array',
                     'set',
                     'set[@comment="ion {}"]',
                     'set[@comment="spin {}"]',
                     "r"])
    path = path.format(self.resort.index(atom_index) + 1, spin)
    log.debug(path)

    results = [[float(x) for x in el.text.split()]
               for el in tree.findall(path)]

    if efermi is None:
        efermi = self.get_fermi_level()
    else:
        efermi = 0.0

    energy = np.array([x[0] for x in results]) - efermi
    ados = np.array([x['spd'.index(orbital) + 1] for x in results])

    return [energy, ados]


@monkeypatch_class(vasp.Vasp)
def get_elapsed_time(self):
    """Return elapsed calculation time in seconds from the OUTCAR file."""
    import re
    regexp = re.compile('Elapsed time \(sec\):\s*(?P<time>[0-9]*\.[0-9]*)')

    with open(os.path.join(self.directory, 'OUTCAR')) as f:
        lines = f.readlines()

    # fragile but fast.
    m = re.search(regexp, lines[-8])

    time = m.groupdict().get('time', None)
    if time is not None:
        return float(time)
    else:
        return None
