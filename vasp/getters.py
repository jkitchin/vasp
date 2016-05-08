import os
import numpy as np
from xml.etree import ElementTree
import vasp
from vasp import log
from monkeypatch import monkeypatch_class

@monkeypatch_class(vasp.Vasp)
def get_k_point_weights(self):
    """Return the k-point weights."""
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
    with open(os.path.join(self.directory,
                           'vasprun.xml')) as f:
        tree = ElementTree.parse(f)
        path = '/'.join(['calculation',
                         'eigenvalues',
                         'array',
                         'set',
                         "set[@comment='spin {}']",
                         "set[@comment='kpoint {}']"])

        # Vasp seems to start at 1 not 0
        fields = tree.find(path.format(spin + 1, kpt + 1))

        return np.array([float(x.text.split()[0]) for x in fields])


@monkeypatch_class(vasp.Vasp)
def get_fermi_level(self):
    """Return the Fermi level."""
    with open(os.path.join(self.directory,
                           'vasprun.xml')) as f:
        tree = ElementTree.parse(f)
        path = '/'.join(['calculation',
                         'dos',
                         "i[@name='efermi']"
                         ])
        return float(tree.find(path).text)
