"""The core Vasp calculator.

I aim to keep this file at a minimum. Hence, many logically grouped
class methods are actually imported at the end.

"""

import os
import subprocess

from ase.calculators.calculator import Calculator
from ase.calculators.calculator import FileIOCalculator

from vasp import log

def VaspExceptionHandler(exc_type, exc_value, exc_traceback):
    """Hande EXC."""
    print('Caught ', exc_type, exc_value, exc_traceback)

class Vasp(FileIOCalculator, object):
    """Class for doing VASP calculations.

    set $ASE_VASP_COMMAND to the command used to run vasp.

    POTCARs are found in:
    $VASP_PP_PATH/potpaw_LDA
    $VASP_PP_PATH/potpaw_PBE
    $VASP_PP_PATH/potpaw_GGA

    """

    name = 'VASP'
    command = None

    implemented_properties = ['energy', 'forces', 'stress',
                              'magmom',  # the overall magnetic moment
                              'magmoms']  # the individual magnetic moments

    # These allow you to use simple strings for the xc kwarg and automatically
    # set the relevant vasp tags.
    xc_defaults = {'lda': {'pp': 'LDA'},
                   # GGAs
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
                   'optb88-vdw': {'pp': 'LDA', 'gga': 'BO',  'luse_vdw': True,
                                  'aggac': 0.0, 'param1': 1.1/6.0,
                                  'param2': 0.22},
                   'optb86b-vdw': {'pp': 'LDA', 'gga': 'MK', 'luse_vdw': True,
                                   'aggac': 0.0, 'param1': 0.1234,
                                   'param2': 1.0},
                   'vdw-df2': {'pp': 'LDA', 'gga': 'ML', 'luse_vdw': True,
                               'aggac': 0.0, 'zab_vdw': -1.8867},
                   'beef-vdw': {'pp': 'LDA', 'gga': 'BF', 'luse_vdw': True,
                                'zab_vdw': -1.8867},
                   # hybrids
                   'hse03': {'pp': 'LDA', 'gga': 'PE', 'lhfcalc': True,
                             'hfscreen': 0.3},
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
        lcharge=False,
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

    def __init__(self, label,
                 restart=None, ignore_bad_restart_file=False,
                 atoms=None, scratch=None,
                 debug=None,
                 exception_handler=VaspExceptionHandler,
                 **kwargs):
        """Create a Vasp calculator.

        label: the directory where the calculation files will be and
        the calculation run.

        debug: an integer, but usually something like logging.DEBUG

        exception_handler: a function that takes single argument of an
        exception, and handles it.

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
        self.exception_handler = exception_handler

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        if debug is not None:
            log.setLevel(debug)

        self.read_metadata()
        log.debug('init metadata = {}'.format(self.metadata))
        # Resort atoms and generate list of POTCARS/counts we need.
        self.sort_atoms(atoms)

        # We need to handle some special kwargs here. The idea is that
        # self.parameters is ready to be written out later. This has to be done
        # after the sorting is known so that arrays of properties can be made
        # the right size and in the right order

        # Spin-polarization.
        if 'ispin' in self.parameters:
            d = self.set_ispin_dict(self.parameters['ispin'])
            self.parameters.update(d)

        # DFT + U uses a dictionary of:
        # {symbol: {'L': int, 'U':float 'J':float}}
        # It is incompatible with setups because the dictionary doesn't allow
        # differentiating between symbols
        if 'ldau_luj' in self.parameters:
            d = self.set_ldau_luj_dict(self.parameters['ldau_luj'])
            self.parameters.update(d)

        # set the exchange-correlation function tags. This sets pp for the
        # potcar path.
        if 'xc' in self.parameters:
            self.parameters['xc'] = self.parameters['xc'].lower()
        self.parameters.update(self.xc_defaults[self.parameters['xc'].lower()])

        # This is opinionated. But necessary for smart
        # reruns. Basically, if a calculation is finished we want to
        # make the atoms consistent with it.
        if os.path.exists(self.outcar):
            with open(self.outcar) as f:
                lines = f.readlines()
                if 'Voluntary context switches:' in lines[-1]:
                    # sets results and updates the atoms.
                    self.read_results()

        if atoms:
            atoms.set_calculator(self)

        # Done with initialization

    def sort_atoms(self, atoms):
        """Generate resort list, and make list of POTCARs to use.

        Returns None.

        """
        self.resort = None
        self.ppp_list = None
        self.symbol_count = None

        if atoms is None:
            return
        self.atoms = atoms
        # Now we sort the atoms and generate the list of POTCARS
        # We end up with ppp = [(index_or_symbol, potcar_file, count)]
        # and resort_indices
        setups = self.parameters.get('setups', [])
        pp = self.parameters['pp']

        ppp = []  # [(index_or_symbol, potcar_file, count)]

        # indices of original atoms needed to make sorted atoms list
        resort_indices = []

        # First the numeric index setups
        for setup in [x for x in setups if isinstance(x[0], int)]:
            ppp += [[setup[0],
                     'potpaw_{}/{}{}/POTCAR'.format(pp, atoms[setup[0]].symbol,
                                                    setup[1]),
                     1]]
            resort_indices += [setup[0]]

        # now the rest of the setups. These are atom symbols
        for setup in [x for x in setups if not isinstance(x[0], int)]:
            symbol = setup[0]
            count = 0
            for i, atom in enumerate(atoms):
                if atom.symbol == symbol and i not in resort_indices:
                    count += 1
                    resort_indices += [i]
            ppp += [[atom.symbol,
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
                if atom.symbol == symbol and i not in resort_indices:
                    resort_indices += [i]
                    count += 1
            if count > 0:
                ppp += [[symbol,
                         'potpaw_{}/{}/POTCAR'.format(pp, symbol),
                         count]]

        assert len(resort_indices) == len(atoms), \
            'Sorting error. sort_indices={}'.format(resort_indices)

        assert sum([x[2] for x in ppp]) == len(atoms)

        self.resort = resort_indices
        self.ppp_list = ppp
        self.atoms_sorted = atoms[self.resort]
        self.symbol_count = [(x[0] if isinstance(x[0], str)
                              else atoms[x[0]].symbol,
                              x[2]) for x in ppp]

        self.metadata['resort'] = self.resort

    def __str__(self):
        """Pretty representation of a calculation.

        TODO: make more like jaspsum.

        """
        s = ['']
        s += ['Vasp calculation in {self.directory}\n']
        if os.path.exists(self.incar):
            with open(self.incar) as f:
                s += [f.read()]
        else:
            s += ['No INCAR yet']

        if os.path.exists(self.poscar):
            with open(self.poscar) as f:
                s += [f.read()]
        else:
            s += ['No POSCAR yet']

        return '\n'.join(s).format(self=self)

    def set_label(self, label):
        """Set working directory.

        In VASP there is no prefix, only the working directory.

        """

        if label is None:
            self.directory = "."
            self.prefix = None
        else:
            self.directory, self.prefix = label, None

        # Convenient attributes for file names
        for f in ['INCAR', 'POSCAR', 'CONTCAR', 'POTCAR',
                  'KPOINTS', 'OUTCAR']:
            fname = os.path.join(self.directory, f)
            setattr(self, f.lower(), fname)

    def check_state(self, atoms):
        """TODO: what should this do?"""
        system_changes = FileIOCalculator.check_state(self, atoms)
        # Ignore boundary conditions:
        if 'pbc' in system_changes:
            system_changes.remove('pbc')

        return system_changes

    def set(self, **kwargs):
        """Set parameters with keyword=value pairs.

        calc.set(xc='PBE')

        A few special kwargs are handled separately to expand them
        prior to setting the parameters. This is done to enable one
        set to track changes.

        """

        if 'xc' in kwargs:
            kwargs.update(self.set_xc_dict(kwargs['xc']))

        if 'ispin' in kwargs:
            kwargs.update(self.set_ispin_dict(kwargs['ispin']))

        if 'ldau_luj' in kwargs:
            kwargs.update(self.set_ldau_luj_dict(kwargs['ldau_luj']))

        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()
        return changed_parameters

    def reset(self):
        """overwrite to avoid killing self.atoms."""
        self.results = {}

    def read_results(self):
        """Read energy, forces from output file.

        Other quantities will be read by other functions.

        """
        from ase.io.vasp import read_vasp_xml

        atoms = read_vasp_xml(os.path.join(self.directory,
                                           'vasprun.xml')).next()

        # update the atoms
        self.atoms = atoms

        self.results['energy'] = atoms.get_potential_energy()
        self.results['forces'] = atoms.get_forces()[self.resort]

    def calculation_required(self, atoms, properties):
        """Returns if a calculation is needed."""

        # if the calculation is finished we do not need to run.
        if os.path.exists(self.outcar):
            with open(self.outcar) as f:
                lines = f.readlines()
                if 'Voluntary context switches:' in lines[-1]:
                    return False

        system_changes = self.check_state(atoms)
        if system_changes:
            return True

        for name in properties:
            if name not in self.results:
                return True

        # default to not required.
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
