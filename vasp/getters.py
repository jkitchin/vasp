import os
import re

import numpy as np
from xml.etree import ElementTree
import vasp
from vasp import log
from monkeypatch import monkeypatch_class


@monkeypatch_class(vasp.Vasp)
def get_ibz_k_points(self):
    lines = open(os.path.join(self.directory, 'OUTCAR'), 'r').readlines()
    ibz_kpts = []
    n = 0
    i = 0
    for line in lines:
        if line.rfind('Following cartesian coordinates') > -1:
            m = n + 2
            while i == 0:
                ibz_kpts.append([float(lines[m].split()[p])
                                 for p in range(3)])
                m += 1
                if lines[m] == ' \n':
                    i = 1
        if i == 1:
            continue
        n += 1
    ibz_kpts = np.array(ibz_kpts)
    return np.array(ibz_kpts)


@monkeypatch_class(vasp.Vasp)
def get_occupation_numbers(self, kpt=0, spin=0):
    lines = open(os.path.join(self.directory, 'OUTCAR')).readlines()
    nspins = self.get_number_of_spins()
    start = 0
    if nspins == 1:
        for n, line in enumerate(lines):  # find it in the last iteration
            m = re.search(' k-point *' + str(kpt + 1) + ' *:', line)
            if m is not None:
                start = n
    else:
        for n, line in enumerate(lines):
            # find it in the last iteration
            if line.find(' spin component ' + str(spin + 1)) != -1:
                start = n
        for n2, line2 in enumerate(lines[start:]):
            m = re.search(' k-point *' + str(kpt + 1) + ' *:', line2)
            if m is not None:
                start = start + n2
                break
    for n2, line2 in enumerate(lines[start + 2:]):
        if not line2.strip():
            break
    occ = []
    for line in lines[start + 2:start + 2 + n2]:
        occ.append(float(line.split()[2]))
    return np.array(occ)

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


@monkeypatch_class(vasp.Vasp)
def get_default_number_of_electrons(self, filename=None):
    """Return the default electrons for each species."""
    if filename is None:
        filename = os.path.join(self.directory, 'POTCAR')

    if not os.path.exists(filename):
        self.write_input()
        
    nelect = []
    lines = open(filename).readlines()
    for n, line in enumerate(lines):
        if line.find('TITEL') != -1:
            symbol = line.split('=')[1].split()[1].split('_')[0].strip()
            valence = float(lines[n + 4].split(';')[1]
                            .split('=')[1].split()[0].strip())
            nelect.append((symbol, valence))
    return nelect


@monkeypatch_class(vasp.Vasp)
def get_valence_electrons(self):
    """Return the number of valence electrons for the atoms.
    Calculated from the POTCAR file.
    """

    default_electrons = self.get_default_number_of_electrons()

    d = {}
    for s, n in default_electrons:
        d[s] = n
    atoms = self.get_atoms()

    nelectrons = 0
    for atom in atoms:
        nelectrons += d[atom.symbol]
    return nelectrons
