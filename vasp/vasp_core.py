"""The core Vasp calculator.

I aim to keep this file at a minimum. Hence, many logically grouped
class methods are actually imported at the end.

"""

import logging
import os
import subprocess
import warnings
import numpy as np
import ase
from ase.calculators.calculator import Calculator
from ase.calculators.calculator import FileIOCalculator
from ase.io import read

# internal modules
import exceptions
import validate
from vasprc import VASPRC
from vasp import log


def VaspExceptionHandler(calc, exc_type, exc_value, exc_traceback):
    """Handle exceptions."""
    if exc_type == exceptions.VaspSubmitted:
        print(exc_value)
        return None
    elif exc_type == exceptions.VaspQueued:
        print(exc_value)
        return None
    elif exc_type == KeyError and exc_value.message == 'energy':
        return None
    elif exc_type == KeyError and exc_value.message == 'forces':
        return np.array([[None, None, None] for atom in calc.get_atoms()])
    elif exc_type == KeyError and exc_value.message == 'stress':
        return np.array([None, None, None, None, None, None])

    print('Unhandled exception in Vasp')
    import traceback
    import sys
    traceback.print_exception(exc_type, exc_value, exc_traceback,
                              file=sys.stdout)
    raise


class Vasp(FileIOCalculator, object):
    """Class for doing VASP calculations.

    Configurations are in vasp.vasprc

    POTCARs are found in:
    $VASP_PP_PATH/potpaw_LDA
    $VASP_PP_PATH/potpaw_PBE
    $VASP_PP_PATH/potpaw_GGA

    """
    version = "0.9.3"
    name = 'VASP'
    command = None
    debug = None

    # List of calculators created
    calculators = []

    implemented_properties = ['energy', 'forces', 'stress',
                              'charges', 'dipole',
                              'magmom',  # the overall magnetic moment
                              'magmoms']  # the individual magnetic moments

    # These allow you to use simple strings for the xc kwarg and automatically
    # set the relevant vasp tags.
    xc_defaults = {'lda': {'pp': 'LDA'},
                   # GGAs
                   'gga': {'pp': 'GGA'},
                   'pbe': {'pp': 'PBE'},
                   'revpbe': {'pp': 'LDA', 'gga': 'RE'},
                   'rpbe': {'pp': 'LDA', 'gga': 'RP'},
                   'am05': {'pp': 'LDA', 'gga': 'AM'},
                   'pbesol': {'pp': 'LDA', 'gga': 'PS'},
                   # Meta-GGAs
                   'tpss': {'pp': 'PBE', 'metagga': 'TPSS'},
                   'revtpss': {'pp': 'PBE', 'metagga': 'RTPSS'},
                   'm06l': {'pp': 'PBE', 'metagga': 'M06L'},
                   # vdW-DFs
                   'optpbe-vdw': {'pp': 'LDA', 'gga': 'OR', 'luse_vdw': True,
                                  'aggac': 0.0},
                   'optb88-vdw': {'pp': 'LDA', 'gga': 'BO', 'luse_vdw': True,
                                  'aggac': 0.0, 'param1': 1.1 / 6.0,
                                  'param2': 0.22},
                   'optb86b-vdw': {'pp': 'LDA', 'gga': 'MK', 'luse_vdw': True,
                                   'aggac': 0.0, 'param1': 0.1234,
                                   'param2': 1.0},
                   'vdw-df2': {'pp': 'LDA', 'gga': 'ML', 'luse_vdw': True,
                               'aggac': 0.0, 'zab_vdw': -1.8867},
                   'beef-vdw': {'pp': 'PBE', 'gga': 'BF', 'luse_vdw': True,
                                'zab_vdw': -1.8867, 'lbeefens': True},
                   # hybrids
                   'pbe0': {'pp': 'LDA', 'gga': 'PE', 'lhfcalc': True},
                   'hse03': {'pp': 'LDA', 'gga': 'PE', 'lhfcalc': True,
                             'hfscreen': 0.3},
                   'hse06': {'pp': 'LDA', 'gga': 'PE', 'lhfcalc': True,
                             'hfscreen': 0.2},
                   'b3lyp': {'pp': 'LDA', 'gga': 'B3', 'lhfcalc': True,
                             'aexx': 0.2, 'aggax': 0.72,
                             'aggac': 0.81, 'aldac': 0.19},
                   'hf': {'pp': 'PBE', 'lhfcalc': True, 'aexx': 1.0,
                          'aldac': 0.0, 'aggac': 0.0}}

    default_parameters = dict(
        xc='PBE',
        pp='PBE',
        ismear=1,
        sigma=0.1,
        lwave=False,
        lcharg=False,
        kpts=[1, 1, 1])

    # These need to be kept separate for writing the incar.
    special_kwargs = ['xc',  # sets vasp tags for the exc-functional
                      'pp',  # determines where POTCARs are retrieved from
                      'setups',
                      # kpoints
                      'kpts',
                      'gamma',
                      'kpts_nintersections',
                      'reciprocal',
                      # DFT + U dictionary
                      'ldau_luj']

    # enumerated states
    EMPTY = 0
    NEW = 1
    QUEUED = 2
    FINISHED = 3
    NOTFINISHED = 4
    EMPTYCONTCAR = 5
    NEB = 10
    UNKNOWN = 100

    VASPRC = VASPRC
    log = log

    def __init__(self, label,
                 restart=True, ignore_bad_restart_file=False,
                 atoms=None, scratch=None,
                 debug=None,
                 exception_handler=VaspExceptionHandler,
                 **kwargs):
        """Create a Vasp calculator.

        label: the directory where the calculation files will be and
        the calculation run.

        debug: an integer, but usually something like logging.DEBUG

        exception_handler: A function for
        handling exceptions. The function should take the arguments
        returned by sys.exc_info(), which is the exception type, value
        and traceback. The default is VaspExceptionHandler.

        **kwargs
          Any Vasp keyword can be used, e.g. encut=450.

          The tag will be upcased when written, and the value is
          written depending on its type. E.g. integers, floats and
          strings are written as they are. True/False is written as
          .TRUE. and .FALSE. and Python lists/tuples are written as
          space delimited lists.

        Special kwargs:

        xc: string indicating the functional to use. It is expanded
        from Vasp.xc_defaults to the relevant Vasp tags.

        kpts: Usually a 3 element list of [k1, k2, k3], but may also
        be a list of kpts.

        setups: This describes special setups for the POTCARS. It is a list of
          the following items.

          (atom_index, suffix)   for exampe: (2, '_sv')

          (atom_symbol, suffix)  for example ('Zr', '_sv')

          If (atom_index, suffix) is used then only that atom index will have a
          POTCAR defined by '{}{}'.format(atoms[atom_index].symbol, suffix)

          If (atom_symbol, suffix) is used then atoms with that symbol (except
          any identified by (atom_index, suffix) will use a POTCAR defined by
          '{}{}'.format(atom_symbol, suffix)

          This syntax has changed from the old dictionary format. The
          reason for this is that this sorting must be
          deterministic. Getting keys in a dictionary is not
          deterministic.

        ldau_luj: This is a dictionary to set the DFT+U tags. For
        example, to put U=4 on the d-orbitals (L=2) of Cu, and nothing
        on the oxygen atoms in a calculation use:

            ldau_luj={'Cu':{'L':2,  'U':4.0, 'J':0.0},
                      'O':{'L':-1, 'U':0.0, 'J':0.0}},

        """
        # set first so self.directory is right
        self.set_label(label)
        self.debug = debug
        if debug is not None:
            log.setLevel(debug)
        self.exception_handler = exception_handler

        self.neb = None
        # We have to check for the type here this because an NEB uses
        # a list of atoms objects. We set pbc to be True because that
        # is what is read in from files, and if we don't the atoms
        # look incompatible.
        if atoms is not None and isinstance(atoms, ase.atoms.Atoms):
            atoms.pbc = [True, True, True]
        elif atoms is not None:
            for a in atoms:
                a.pbc = [True, True, True]
            self.neb = True

        if self.neb is not None:
            self.neb = atoms

        if atoms is not None and self.neb is None:
            self.atoms = atoms

        # We do not pass kwargs here. Some of the special kwargs
        # cannot be set at this point since they need to know about
        # the atoms and parameters. This reads params and results from
        # existing files if they are there. It calls self.read().

        if self.neb:
            FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                      str(label))
        else:
            FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                      str(label), atoms)

        # The calculator should be up to date with the file
        # system here.

        # Add default parameters if they aren't set otherwise.
        for key, val in Vasp.default_parameters.iteritems():
            if key not in kwargs and key not in self.parameters:
                kwargs[key] = val

        # Next we update kwargs with the special kwarg
        # dictionaries. ispin, rwigs are special, and needs sorted
        # atoms. so we save it for later.
        if 'ispin' in kwargs:
            ispin = kwargs['ispin']
            del kwargs['ispin']
        else:
            ispin = None

        if 'rwigs' in kwargs:
            rwigs = kwargs['rwigs']
            del kwargs['rwigs']
        else:
            rwigs = None

        if 'ldau_luj' in kwargs:
            ldau_luj = kwargs['ldau_luj']
            del kwargs['ldau_luj']
        else:
            ldau_luj = None

        # Now update the parameters. If there are any new kwargs here,
        # it will reset the calculator and cause a calculation to be
        # run if needed.
        self.set(**kwargs)

        # In case no atoms was on file, and one is passed in, we set
        # it here.

        if self.atoms is None and atoms is not None and self.neb is None:
            self.sort_atoms(atoms)
        elif self.neb is not None:
            self.sort_atoms(self.neb[0])

        # I don't know why this is necessary. but it seems like
        # these get lost and it causes restart issues
        if atoms is not None and self.neb is None:
            aimm = atoms.get_initial_magnetic_moments()
            self.atoms.set_initial_magnetic_moments(aimm)

        # These depend on having atoms already.
        if ispin is not None:
            self.set(**self.set_ispin_dict(ispin))

        if rwigs is not None:
            self.set(**self.set_rwigs_dict(rwigs))

        if ldau_luj is not None:
            self.set(**self.set_ldau_luj_dict(ldau_luj))

        # Finally run validate functions
        if VASPRC['validate']:
            for key, val in self.parameters.iteritems():
                if key in validate.__dict__:
                    f = validate.__dict__[key]
                    f(self, val)
                else:
                    warnings.warn('No validation for {}'.format(key))

        # Store instance in class for future reference
        Vasp.calculators += [self]

        # We define some classmethods that work on the class
        # Here we redefine some instance methods.
        self.run = self._run
        self.abort = self._abort
        self.wait = self._wait

    def sort_atoms(self, atoms=None):
        """Generate resort list, and make list of POTCARs to use.

        Returns None.

        """
        self.resort = None
        self.ppp_list = None
        self.symbol_count = None

        if atoms is None:
            log.debug('Atoms was none.')
            return
        self.atoms = atoms

        # Now we sort the atoms and generate the list of POTCARS
        # We end up with ppp = [(index_or_symbol, potcar_file, count)]
        # and resort_indices
        setups = self.parameters.get('setups', [])
        pp = self.parameters['pp']

        ppp = []  # [(index_or_symbol, potcar_file, count)]

        # indices of original atoms needed to make sorted atoms list
        sort_indices = []

        # First the numeric index setups
        for setup in [x for x in setups if isinstance(x[0], int)]:
            ppp += [[setup[0],
                     'potpaw_{}/{}{}/POTCAR'.format(pp, atoms[setup[0]].symbol,
                                                    setup[1]),
                     1]]
            sort_indices += [setup[0]]

        # now the rest of the setups. These are atom symbols
        for setup in [x for x in setups if not isinstance(x[0], int)]:
            symbol = setup[0]
            count = 0
            for i, atom in enumerate(atoms):
                if atom.symbol == symbol and i not in sort_indices:
                    count += 1
                    sort_indices += [i]

            ppp += [[symbol,
                     'potpaw_{}/{}{}/POTCAR'.format(pp, symbol, setup[1]),
                     count]]
        # now the remaining atoms use default potentials
        # First get the chemical symbols that remain
        symbols = []
        for atom in atoms or []:
            if (atom.symbol not in symbols and
                atom.symbol not in [x[0] for x in ppp]):
                symbols += [atom.symbol]

        for symbol in symbols:
            count = 0
            for i, atom in enumerate(atoms):
                if atom.symbol == symbol and i not in sort_indices:
                    sort_indices += [i]
                    count += 1
            if count > 0:
                ppp += [[symbol,
                         'potpaw_{}/{}/POTCAR'.format(pp, symbol),
                         count]]

        assert len(sort_indices) == len(atoms), \
            'Sorting error. sort_indices={}'.format(sort_indices)

        assert sum([x[2] for x in ppp]) == len(atoms)
        self.sort = sort_indices

        # This list is used to convert Vasp ordering back to the
        # user-defined order.
        self.resort = [k[1] for k in
                       sorted([[j, i] for i, j in enumerate(sort_indices)])]

        # June 23, 2016. Jake Boes found a bug in how sorting is
        # done. We fixed it, but the fix is not backwards compatible
        # with the old resort we stored in DB.db. It appears we can
        # detect an old case if there is inconsistency with what is
        # stored and calculated here.  We check here to see if we are
        # consistent, and if not fix the issue. This should only occur
        # once.
        if (self.resort is not None
            and self.get_db('resort') is not None
            and self.resort != list(self.get_db('resort'))):
            ns =  [k[1] for k in
                   sorted([[j, i]
                           for i, j in enumerate(self.get_db('resort'))])]
            from ase.db import connect
            with connect(os.path.join(self.directory, 'DB.db')) as con:
                tatoms = con.get_atoms(id=1)
            self.write_db(atoms=tatoms, data={'resort': ns})
            print('Fixed resort issue in {}. '
                  'You should not see this message'
                  ' again'.format(self.directory))
            self.resort = ns
            sort_indices = [k[1] for k in
                            sorted([[j, i]
                                    for i, j in enumerate(ns)])]

        self.ppp_list = ppp
        self.atoms_sorted = atoms[sort_indices]
        self.symbol_count = [(x[0] if isinstance(x[0], str)
                              else atoms[x[0]].symbol,
                              x[2]) for x in ppp]


    def _repr_html_(self):
        """Output function for Jupyter notebooks."""
        from ase.io import write
        atoms = self.get_atoms()
        atoms_image = os.path.join(self.directory, '_repr_html.png')
        path = os.path.relpath(atoms_image, os.getcwd())
        write(atoms_image,
              atoms, show_unit_cell=2)
        formula = atoms.get_chemical_formula()
        energy = self.results.get('energy', None)
        template = """
        <table><tr>
        <td><img src="{path}"></td>
        <td>
        Formula = {formula} <br>
        Energy = {energy} eV   </td>
        </tr>
        </table>
        """

        return template.format(**locals())

    def __str__(self):
        """Pretty representation of a calculation.

        TODO: 1. Incorporate magmoms and vibrational frequencies.
              2. Implement convergence check?

        """
        s = ['\n', 'Vasp calculation directory:']
        s += ['---------------------------']
        s += ['  [[{self.directory}]]']

        atoms = self.get_atoms()
        cell = atoms.get_cell()

        A, B, C = [i for i in cell]
        l = map(np.linalg.norm, cell)
        a, b, c = l
        alpha = np.arccos(np.dot(B / b, C / c)) * 180 / np.pi
        beta = np.arccos(np.dot(A / a, C / c)) * 180 / np.pi
        gamma = np.arccos(np.dot(A / a, B / b)) * 180 / np.pi

        # Format unit cell output
        #########################
        s += ['\nUnit cell:']
        s += ['----------']
        s += ['    {:^8}{:^8}{:^8}{:>12}'.format('x', 'y', 'z',
                                                 '|v|')]
        for i, v in enumerate(cell):
            s += ['  v{0}{2:>8.3f}{3:>8.3f}{4:>8.3f}'
                  '{1:>12.3f} Ang'.format(i, l[i], *v)]

        s += ['  alpha, beta, gamma (deg):'
              '{:>6.1f}{:>6.1f}{:>6.1f}'.format(alpha, beta, gamma)]

        volume = atoms.get_volume()
        s += ['  Total volume:{:>25.3f} Ang^3'.format(volume)]

        # Format stress output
        #########################
        stress = self.results.get('stress', [np.nan] * 6)

        s += ['  Stress:{:>6}{:>7}{:>7}'
              '{:>7}{:>7}{:>7}'.format('xx', 'yy', 'zz',
                                       'yz', 'xz', 'xy')]
        s += ['{:>15.3f}{:7.3f}{:7.3f}'
              '{:7.3f}{:7.3f}{:7.3f} GPa\n'.format(*stress)]

        # Format atoms output
        #########################
        s += ['  {:<4}{:<8}{:<3}{:^10}{:^10}{:^10}'
              '{:>14}'.format('ID', 'tag', 'sym',
                              'x', 'y', 'z', 'rmsF (eV/A)')]

        from ase.constraints import FixAtoms, FixScaled
        constraints = [[None, None, None] for atom in atoms]
        for constraint in atoms.constraints:
            if isinstance(constraint, FixAtoms):
                for i in constraint.index:
                    constraints[i] = [True, True, True]
            elif isinstance(constraint, FixScaled):
                constraints[constraint.a] = constraint.mask.tolist()

        forces = self.results.get('forces', np.array([[np.nan, np.nan, np.nan]
                                                      for atom in self.atoms]))
        for i, atom in enumerate(atoms):
            rms_f = np.sum(forces[i] ** 2) ** 0.5

            s += ['  {:<4}{:<8}{:3}{:9.3f}{}{:9.3f}{}{:9.3f}{}'
                  '{:>10.2f}'.format(i, atom.tag, atom.symbol,
                                     atom.x,
                                     '*' if constraints[i][0]
                                     is True else ' ',
                                     atom.y,
                                     '*' if constraints[i][1]
                                     is True else ' ',
                                     atom.z,
                                     '*' if constraints[i][1]
                                     is True else ' ',
                                     rms_f)]

        energy = self.results.get('energy', np.nan)
        s += ['  Potential energy: {:.4f} eV'.format(energy)]

        # Format INPUT output
        #########################
        s += ['\nINPUT Parameters:']
        s += ['-----------------']
        for key, value in self.parameters.iteritems():
            s += ['  {0:<10}: {1}'.format(key, value)]

        # Format pseudo-potential output
        #########################
        s += ['\nPseudopotentials used:']
        s += ['----------------------']
        for sym, ppp, hash in self.get_pseudopotentials():
            s += ['  {}: {} (git-hash: {})'.format(sym, ppp, hash)]

        return '\n'.join(s).format(self=self)

    def set_label(self, label):
        """Set working directory.

        In VASP there is no prefix, only the working directory.

        """

        if label is None:
            self.directory = os.path.abspath(".")
            self.prefix = None
        else:
            d = os.path.expanduser(label)
            d = os.path.abspath(d)
            self.directory, self.prefix = d, None
            if not os.path.isdir(self.directory):
                os.makedirs(self.directory)

        # Convenient attributes for file names
        for f in ['INCAR', 'POSCAR', 'CONTCAR', 'POTCAR',
                  'KPOINTS', 'OUTCAR']:
            fname = os.path.join(self.directory, f)
            setattr(self, f.lower(), fname)

    def check_state(self, atoms=None):
        """Check if any changes exist that require new calculations."""
        if atoms is None:
            atoms = self.get_atoms()

        log.debug('atoms IMM: {}'.format(atoms.get_initial_magnetic_moments()))
        system_changes = FileIOCalculator.check_state(self, atoms)
        # Ignore boundary conditions:
        if 'pbc' in system_changes:
            system_changes.remove('pbc')

        s = 'FileIOCalculator reports these changes: {}'
        log.debug(s.format(system_changes))
        # if dir is empty, there is nothing to read here.
        if self.get_state() == Vasp.EMPTY:
            return system_changes

        # Check if the parameters have changed
        file_params = {}
        file_params.update(self.read_incar())
        file_params.update(self.read_potcar())
        file_params.update(self.read_kpoints())

        xc_keys = sorted(Vasp.xc_defaults,
                         key=lambda k: len(Vasp.xc_defaults[k]),
                         reverse=True)

        for ex in xc_keys:
            pd = {k: file_params.get(k, None)
                  for k in Vasp.xc_defaults[ex]}
            if pd == Vasp.xc_defaults[ex]:
                file_params['xc'] = ex.lower()
                break

        # reconstruct ldau_luj if necessary
        if 'ldauu' in file_params:
            ldaul = file_params['ldaul']
            ldauj = file_params['ldauj']
            ldauu = file_params['ldauu']

            with open(self.potcar) as f:
                lines = f.readlines()

            # symbols are in the first line of each potcar
            symbols = [lines[0].split()[1]]
            for i, line in enumerate(lines):
                if 'End of Dataset' in line and i != len(lines) - 1:
                    symbols += [lines[i + 1].split()[1]]

            ldau_luj = {}
            for sym, l, j, u in zip(symbols, ldaul, ldauj, ldauu):
                ldau_luj[sym] = {'L': l, 'U': u, 'J': j}

            file_params['ldau_luj'] = ldau_luj

        if not {k: v for k, v in self.parameters.iteritems()
                if v is not None} == file_params:
            new_keys = set(self.parameters.keys()) - set(file_params.keys())
            missing_keys = (set(file_params.keys()) -
                            set(self.parameters.keys()))
            log.debug('New keys: {}'.format(new_keys))
            log.debug('Missing keys: {}'.format(missing_keys))
            log.debug('params_on_file do not match.')
            log.debug('file-params: {}'.format(file_params))
            log.debug('compared to: {}'.format({k: v for k, v in
                                                self.parameters.iteritems()
                                                if v is not None}))
            system_changes += ['params_on_file']

        log.debug('System changes: {}'.format(system_changes))
        return system_changes

    def reset(self):
        """overwrite to avoid killing self.atoms."""
        log.debug('Resetting calculator.')
        self.results = {}

    def update(self, atoms=None):
        """Updates calculator.

        If a calculation is required,  run it, otherwise updates results.

        """
        if atoms is None:
            atoms = self.get_atoms()

        if self.neb:
            return self.get_neb()

        if self.calculation_required(atoms, ['energy']):
            return self.calculate(atoms)
        else:
            self.read_results()

        return True

    def calculation_required(self, atoms=None, properties=['energy']):
        """Returns if a calculation is needed."""

        if atoms is None:
            atoms = self.get_atoms()

        system_changes = self.check_state(atoms)
        if system_changes:
            log.debug('Calculation needed for {}'.format(system_changes))
            return True
        for name in properties:
            if name not in self.results:
                log.debug('{} not in {}. Calc required.'.format(name,
                                                                self.results))
                return True

        # if the calculation is finished we do not need to run.
        if os.path.exists(self.outcar):
            with open(self.outcar) as f:
                lines = f.readlines()
                if 'Voluntary context switches:' in lines[-1]:
                    return False

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=None):
        """Runs a calculation, only if necessary."""
        if self.calculation_required(atoms, properties):

            # The subclass implementation should first call this
            # implementation to set the atoms attribute.
            Calculator.calculate(self, atoms, properties, system_changes)

            self.write_input(atoms, properties, system_changes)

            if self.command is None:
                raise RuntimeError('Please set $%s environment variable ' %
                                   ('ASE_' + self.name.upper() + '_COMMAND') +
                                   'or supply the command keyword')

            olddir = os.getcwd()
            try:
                os.chdir(self.directory)
                errorcode = subprocess.call(self.command,
                                            stdout=subprocess.PIPE,
                                            shell=True)

            finally:
                os.chdir(olddir)

            if errorcode:
                s = '{} returned an error: {}'
                raise RuntimeError(s.format(self.name, errorcode))

        # This sets self.results, and updates the atoms
        self.read_results()

    def clone(self, newdir):
        """Copy the calculation directory to newdir and set label to
        newdir.

        """
        state = self.get_state()

        import shutil
        if not os.path.isdir(newdir):
            shutil.copytree(self.directory, newdir)

            # need some cleanup here. do not copy jobids, etc...
            # What survives depends on the state
            # delete these files if not finished.
            if state in [Vasp.QUEUED, Vasp.NOTFINISHED]:
                os.unlink(os.path.join(newdir, 'OUTCAR'))
                os.unlink(os.path.join(newdir, 'vasprun.xml'))

            if state in [Vasp.EMPTYCONTCAR]:
                os.unlink(os.path.join(newdir, 'OUTCAR'))
                os.unlink(os.path.join(newdir, 'vasprun.xml'))
                os.unlink(os.path.join(newdir, 'CONTCAR'))

            # eliminate jobid on copying.
            newdb = os.path.join(newdir, 'DB.db')
            self.write_db(fname=newdb,
                          data={'path': os.path.abspath(newdir),
                                'jobid': None})

        self.__init__(newdir)

    def get_state(self):
        """Determine calculation state based on directory contents.

        Returns an integer for the state.

        """
        # We do not check for KPOINTS here. That file may not exist if
        # the kspacing incar parameter is used.
        base_input = [os.path.exists(os.path.join(self.directory, f))
                       for f in ['INCAR', 'POSCAR', 'POTCAR']]

        # Check for NEB first.
        if (np.array([os.path.exists(os.path.join(self.directory, f))
                      for f in ['INCAR', 'POTCAR']]).all()
            and not os.path.exists(os.path.join(self.directory, 'POSCAR'))
            and os.path.isdir(os.path.join(self.directory, '00'))):
            return Vasp.NEB

        # Some input does not exist
        if False in base_input:
            # some input file is missing
            return Vasp.EMPTY

        # Input files exist, but no jobid, and no output
        if (np.array(base_input).all()
            and self.get_db('jobid') is not None
            and not os.path.exists(os.path.join(self.directory, 'OUTCAR'))):
            return Vasp.NEW

        # INPUT files exist, a jobid in the queue
        if self.in_queue():
            return Vasp.QUEUED

        # Not in queue, and finished
        if not self.in_queue():
            if os.path.exists(self.outcar):
                with open(self.outcar) as f:
                    lines = f.readlines()
                    if 'Voluntary context switches:' in lines[-1]:
                        return Vasp.FINISHED

        # Not in queue, and not finished
        if not self.in_queue():
            if os.path.exists(self.outcar):
                with open(self.outcar) as f:
                    lines = f.readlines()
                    if 'Voluntary context switches:' not in lines[-1]:
                        return Vasp.NOTFINISHED
            else:
                return Vasp.NOTFINISHED

        # Not in queue, and not finished, with empty contcar
        if not self.in_queue():
            if os.path.exists(self.contcar):
                with open(self.contcar) as f:
                    if f.read() == '':
                        return Vasp.EMPTYCONTCAR

        return Vasp.UNKNOWN

    @property
    def potential_energy(self):
        """Property to return potential_energy."""
        self.update()
        atoms = self.get_atoms()
        return atoms.get_potential_energy()

    @property
    def forces(self, apply_constraints=False):
        """Property to return forces."""
        self.update()
        atoms = self.get_atoms()
        return atoms.get_forces(apply_constraints)

    @property
    def stress(self):
        """Property to return stress."""
        self.update()
        atoms = self.get_atoms()
        return atoms.get_stress()

    @property
    def traj(self):
        """Get a trajectory.

        This reads Atoms objects from vasprun.xml. By default returns
        all images.  If index is an integer, return that image.

        Technically, this is just a list of atoms with a
        SinglePointCalculator attached to them.

        This is usually only relevant if you have done a
        relaxation. If the calculation is an NEB, the images are
        returned.

        """
        from ase.calculators.singlepoint import SinglePointCalculator as SPC
        self.update()

        if self.neb:
            images, energies = self.get_neb()
            tatoms = [x.copy() for x in images]
            for i, x in enumerate(tatoms):
                x.set_calculator(SPC(x, energy=energies[i]))
            return tatoms

        LOA = []
        for atoms in read(os.path.join(self.directory, 'vasprun.xml'), ':'):
            catoms = atoms.copy()
            catoms = catoms[self.resort]
            catoms.set_calculator(SPC(catoms,
                                      energy=atoms.get_potential_energy()))
            LOA += [catoms]
        return LOA

    def view(self, index=None):
        """Visualize the calculation.

        """
        from ase.visualize import view
        if index is not None:
            return view(self.traj[index])
        else:
            return view(self.traj)

    def describe(self, long=False):
        """Describe the parameters used with docstrings in vasp.validate."""
        for key in sorted(self.parameters.keys()):
            if key in validate.__dict__:
                f = validate.__dict__[key]
                d = f.__doc__ or 'No docstring found.'
                print('{} = {}:'.format(key, self.parameters[key]))
                if long:
                    print('  ' + d)
                else:
                    print('  ' + d.split('\n')[0])
                print('')

    @property
    def ready(self):
        """Property for is calculator ready.

        That means no calculation is required to get results. Checking
        this should not trigger a calculation.

        """
        return not self.calculation_required()

    @classmethod
    def run(cls, wait=False):
        """Convenience function to run calculators.

        The default behavior is to exit after doing this. If wait is
        True, iy will cause it to wait with the default args to
        Vasp.wait.

        If wait is a dictionary, it will be passed as kwargs to
        Vasp.wait.

        """
        energies = [calc.potential_energy for calc in Vasp.calculators]

        if None not in energies:
            # They are all done.
            return energies

        if wait is False:
            Vasp.abort()
        elif isinstance(wait, dict):
            Vasp.wait(**wait)
        else:
            Vasp.wait()

    def _run(self):
        """Convenience function to run the calculator."""
        return self.potential_energy

    @classmethod
    def all(cls):
        """Returns if all calculators in the class are ready."""
        status = [c.ready for c in Vasp.calculators]
        return False not in status

    @classmethod
    def stop_if(cls, condition):
        """Stops the program if condition is truthy."""
        if condition:
            import sys
            sys.exit()

    @classmethod
    def abort(cls):
        """Abort and exit the program the calculator is running in."""
        Vasp.stop_if(True)

    def _abort(self):
        """Abort and exit program the calculator is running in."""
        Vasp.stop_if(True)

    @classmethod
    def clear_calculators(cls):
        """Clear the stored calculators."""
        Vasp.calculators = []

    @classmethod
    def vasprc(cls, **kwargs):
        """Convenience method to update VASPRC.

        Vasp.vasprc(mode=None)
        """
        Vasp.VASPRC.update(kwargs)

    @classmethod
    def wait(cls, poll_interval=5, timeout=None, abort=False):
        """Control function to wait until all calculators are ready.

        if abort is truthy, stop the program.

        Otherwise check the calculators every poll_interval seconds,
        up to timeout seconds later. If timeout is None, poll forever.

        """

        if abort and not Vasp.all():
            Vasp.abort()

        import time
        if timeout is not None:
            t0 = time.time()
            while not Vasp.all() and time.time() - t0 < timeout:
                time.sleep(poll_interval)

            if time.time() - t0 > timeout:
                print('Timeout exceeded without finishing.')
                Vasp.abort()

        else:
            while not Vasp.all():
                time.sleep(poll_interval)

    def _wait(self, poll_interval=5, timeout=None, abort=False):
        """Control function to wait until all calculators are ready.

        if abort is truthy, stop the program.

        Otherwise check the calculators every poll_interval seconds,
        up to timeout seconds later. If timeout is None, poll forever.

        """
        self.update()
        if abort and not self.ready:
            self.abort()

        import time
        if timeout is not None:
            t0 = time.time()
            while not self.ready and time.time() - t0 < timeout:
                time.sleep(poll_interval)

            if time.time() - t0 > timeout:
                print('Timeout exceeded without finishing.')
                self.abort()

        else:
            while not self.ready:
                time.sleep(poll_interval)
