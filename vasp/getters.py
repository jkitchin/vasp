import os
import re
from hashlib import sha1

import numpy as np
from xml.etree import ElementTree
import vasp
from vasp import log
from monkeypatch import monkeypatch_class


@monkeypatch_class(vasp.Vasp)
def get_beefens(self, n=-1):
    """Get the BEEFens 2000 ensemble energies from the OUTCAR.
    This only works with Vasp 5.3.5 compiled with libbeef.
    I am pretty sure this array is the deviations from the total
    energy. There are usually 2000 of these, but it is not clear this will
    always be the case. I assume the 2000 entries are always in the same
    order, so you can calculate ensemble energy differences for reactions,
    as long as the number of samples in the ensemble is the same.
    There is usually more than one BEEFens section. By default we
    return the last one. Choose another one with the the :par: n.
    see http://suncat.slac.stanford.edu/facility/software/functional/
    """
    self.update()
    beefens = []
    with open(os.path.join(self.directory, 'OUTCAR')) as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if 'BEEFens' in line:
                nsamples = int(re.search('(\d+)', line).groups()[0])
                beefens.append([float(x) for x in lines[i + 1: i + nsamples]])
    return np.array(beefens[n])


@monkeypatch_class(vasp.Vasp)
def get_ibz_k_points(self):
    """Return the irreducible k-points."""
    self.update()
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
    """Return the occupation of each k-point."""
    self.update()
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
    self.update()

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
    self.update()
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
        path = path.format(spin + 1, kpt + 1)
        log.debug('path={}'.format(path))
        # Vasp seems to start at 1 not 0
        fields = tree.find(path)

        return np.array([float(x.text.split()[0]) for x in fields])


@monkeypatch_class(vasp.Vasp)
def get_fermi_level(self):
    """Return the Fermi level."""
    self.update()

    with open(os.path.join(self.directory,
                           'vasprun.xml')) as f:
        tree = ElementTree.parse(f)
        path = '/'.join(['calculation',
                         'dos',
                         "i[@name='efermi']"
                         ])
        return float(tree.find(path).text)


@monkeypatch_class(vasp.Vasp)
def get_ados(self, atom_index, orbital, spin=1, efermi=None):
    """Return Atom projected DOS for atom index, orbital and spin.

    orbital: string ['s', 'p', 'd']

    If efermi is not None, use this value as 0.0.

    :returns: (energies, ados)

    """
    self.update()

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


@monkeypatch_class(vasp.Vasp)
def get_volumetric_data(self, filename=None, **kwargs):
    """Read filename to read the volumetric data in it.
    Supported filenames are CHG, CHGCAR, and LOCPOT.
    """
    if filename is None:
        filename = os.path.join(self.directory, 'CHG')

    from VaspChargeDensity import VaspChargeDensity

    atoms = self.get_atoms()
    vd = VaspChargeDensity(filename)

    data = np.array(vd.chg)
    n0, n1, n2 = data[0].shape

    # This is the old code, but it doesn't seem to work anymore.
    #s0 = np.linspace(0, 1, num=n0, endpoint=False)
    #s1 = np.linspace(0, 1, num=n1, endpoint=False)
    #s2 = np.linspace(0, 1, num=n2, endpoint=False)

    # X, Y, Z = np.meshgrid(s0, s1, s2)

    s0 = 1.0 / n0
    s1 = 1.0 / n1
    s2 = 1.0 / n2
    X, Y, Z = np.mgrid[0.0:1.0:s0,
                       0.0:1.0:s1,
                       0.0:1.0:s2]

    C = np.column_stack([X.ravel(),
                         Y.ravel(),
                         Z.ravel()])

    uc = atoms.get_cell()
    real = np.dot(C, uc)

    # now convert arrays back to unitcell shape
    x = np.reshape(real[:, 0], (n0, n1, n2))
    y = np.reshape(real[:, 1], (n0, n1, n2))
    z = np.reshape(real[:, 2], (n0, n1, n2))
    return (x, y, z, data)


@monkeypatch_class(vasp.Vasp)
def get_charge_density(self, spin=0, filename=None):
    """Returns x, y, and z coordinate and charge density arrays.

    Supported file formats: CHG, CHGCAR
    :param int spin: an integer
    :returns: x, y, z, charge density arrays
    :rtype: 3-d numpy arrays
    Relies on :func:`ase.calculators.vasp.VaspChargeDensity`.
    """
    self.update()

    if not self.parameters.get('lcharg', False):
        raise Exception('CHG was not written. Set lcharg=True')

    if filename is None:
        filename = os.path.join(self.directory, 'CHG')

    x, y, z, data = get_volumetric_data(self, filename=filename)
    return x, y, z, data[spin]


@monkeypatch_class(vasp.Vasp)
def get_local_potential(self):
    """Returns x, y, z, and local potential arrays

    We multiply the data by the volume because we are reusing the
    charge density code which divides by volume.
    """
    self.update()

    fname = os.path.join(self.directory, 'LOCPOT')
    x, y, z, data = get_volumetric_data(self, filename=fname)
    atoms = self.get_atoms()
    return x, y, z, data[0] * atoms.get_volume()


@monkeypatch_class(vasp.Vasp)
def get_elf(self):
    """Returns x, y, z and electron localization function arrays."""

    fname = os.path.join(self.directory, 'ELFCAR')
    x, y, z, data = get_volumetric_data(self, filename=fname)
    atoms = self.get_atoms()
    return x, y, z, data[0] * atoms.get_volume()


@monkeypatch_class(vasp.Vasp)
def get_electron_density_center(self, spin=0, scaled=True):
    """Returns center of electron density.
    If scaled, use scaled coordinates, otherwise use cartesian
    coordinates.
    """

    atoms = self.get_atoms()

    x, y, z, cd = self.get_charge_density(spin)
    n0, n1, n2 = cd.shape
    nelements = n0 * n1 * n2
    voxel_volume = atoms.get_volume() / nelements
    total_electron_charge = cd.sum() * voxel_volume

    electron_density_center = np.array([(cd * x).sum(),
                                        (cd * y).sum(),
                                        (cd * z).sum()])
    electron_density_center *= voxel_volume
    electron_density_center /= total_electron_charge

    if scaled:
        uc = atoms.get_cell()
        return np.dot(np.linalg.inv(uc.T), electron_density_center.T).T
    else:
        return electron_density_center


@monkeypatch_class(vasp.Vasp)
def get_dipole_vector(self, atoms=None):
    """Tries to return the dipole vector of the unit cell in atomic units.

    Returns None when CHG file is empty/not-present.

    """
    self.update()

    from POTCAR import get_ZVAL

    if atoms is None:
        atoms = self.get_atoms()

    try:
        x, y, z, cd = self.get_charge_density()
    except (IOError, IndexError):
        # IOError: no CHG file, function called outside context manager
        # IndexError: Empty CHG file, Vasp run with lcharg=False
        return None

    n0, n1, n2 = cd.shape

    nelements = n0 * n1 * n2
    voxel_volume = atoms.get_volume() / nelements
    total_electron_charge = -cd.sum() * voxel_volume

    electron_density_center = np.array([(cd * x).sum(),
                                        (cd * y).sum(),
                                        (cd * z).sum()])
    electron_density_center *= voxel_volume
    electron_density_center /= total_electron_charge

    electron_dipole_moment = electron_density_center * total_electron_charge
    electron_dipole_moment *= -1.0

    # now the ion charge center
    LOP = self.get_pseudopotentials()
    ppp = os.environ['VASP_PP_PATH']

    # make dictionary for ease of use
    zval = {}
    for sym, ppath, hash in LOP:
        fullpath = os.path.join(ppp,  ppath)
        z = get_ZVAL(fullpath)
        zval[sym] = z

    ion_charge_center = np.array([0.0, 0.0, 0.0])
    total_ion_charge = 0.0
    for atom in atoms:
        Z = zval[atom.symbol]
        total_ion_charge += Z
        pos = atom.position
        ion_charge_center += Z*pos

    ion_charge_center /= total_ion_charge
    ion_dipole_moment = ion_charge_center * total_ion_charge

    dipole_vector = (ion_dipole_moment + electron_dipole_moment)

    return dipole_vector


@monkeypatch_class(vasp.Vasp)
def get_dipole_moment(self, atoms=None):
    """Return dipole_moment.

    dipole_moment = ((dipole_vector**2).sum())**0.5/Debye

    """
    self.update()

    dv = self.get_dipole_vector(atoms)

    from ase.units import Debye
    return ((dv**2).sum())**0.5 / Debye



@monkeypatch_class(vasp.Vasp)
def get_pseudopotentials(self):
    """Return list of (symbol, path, git-hash) for each POTCAR."""
    symbols = [x[0] for x in self.ppp_list]
    paths = [x[1] for x in self.ppp_list]
    hashes = []
    vasp_pp_path = os.environ['VASP_PP_PATH']
    for ppp in paths:
        with open(os.path.join(vasp_pp_path, ppp), 'r') as f:
            data = f.read()

        s = sha1()
        s.update("blob %u\0" % len(data))
        s.update(data)
        hashes.append(s.hexdigest())

    return zip(symbols, paths, hashes)
