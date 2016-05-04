import os
import subprocess

from ase.calculators.calculator import Calculator,\
    FileIOCalculator


class Vasp(FileIOCalculator):
    """Class for doing VASP calculations.

    set $ASE_VASP_COMMAND to the command used to run vasp.

    POTCARs are found in:
    $VASP_PP_PATH/potpaw_LDA
    $VASP_PP_PATH/potpaw_PBE
    $VASP_PP_PATH/potpaw_GGA

    """
    name = 'VASP'
    command = None

    implemented_properties = ['energy', 'forces', 'stress', 'magmom']

    default_parameters = dict(
        xc='PBE',
        ismear=1,
        sigma=0.1,
        lwave=False,
        lcharge=False,
        kpts=[1, 1, 1])

    # These allow you to use simple strings for the xc kwarg and automatically
    # set the relevant vasp tags.
    xc_defaults = {None: {'pp': 'PBE'},
                   'lda': {'pp': 'LDA'},
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

    # These need to be kept separate for writing the incar.
    special_kwargs = ['xc',  # sets vasp tags for the exchange correlation functional
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
                 atoms=None, scratch=None, **kwargs):
        """Create a Vasp calculator.


        setups: This describes special setups for the POTCARS. It is a list of
          the following items.

          (atom_index, suffix)   for exampe: (2, '_sv')

          (atom_symbol, suffix)  for example ('Zr', '_sv')


          If (atom_index, suffix) is used then only that atom index will have a
          POTCAR defined by '{}{}'.format(atoms[atom_index].symbol, suffix)

          If (atom_symbol suffix) is used then atoms with that symbol (except
          any identified by (atom_index, suffix) will use a POTCAR defined by
          '{}{}'.format(atom_symbol, suffix)

          This syntax has changed from the old dictionary format. The reason for
          this is that this sorting must be deterministic. Getting keys in a
          dictionary is not deterministic.

        """

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        # set the exchange-correlation function tags. This sets pp for the
        # potcar path.
        self.parameters.update(self.xc_defaults[self.parameters['xc'].lower()])

        # Now we sort the atoms and generate the list of POTCARS
        # We end up with ppp = [(index_or_symbol, potcar_file, count)]
        # and sort_indices
        setups = self.parameters.get('setups', [])
        pp = self.parameters['pp']

        ppp = []      # [(index_or_symbol, potcar_file, count)]

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
            ppp += [[atom.symbol,
                     'potpaw_{}/{}{}/POTCAR'.format(pp, symbol, setup[1]),
                     count]]
        # now the remaining atoms use default potentials
        # First get the chemical symbols that remain
        symbols = []
        for atom in atoms:
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

        self.sort_indices = sort_indices
        self.ppp_list = ppp
        self.atoms = atoms
        self.atoms_sorted = atoms[sort_indices]
        self.symbol_count = [(x[0] if isinstance(x[0], str)
                              else atoms[x[0]].symbol,
                              x[2]) for x in ppp]

        # TODO: do we need to save this sort data? it used to be in ase-sort.dat
        # if so, it should be saved by a writer.

        # We need to handle some special kwargs here. The idea is that
        # self.parameters is ready to be written out later. This has to be done
        # after the sorting is known so that arrays of properties can be made
        # the right size and in the right order

        # Spin-polarization. there are two ways to get magmoms in.
        # 1. if you use magmoms as a keyword, they are used.
        # 2. if you set magmom on each atom in an Atoms object and do not use
        # magmoms then we use the atoms magmoms, if we have ispin=2 set.
        if self.parameters.get('ispin', None)==2 and not 'magmoms' in self.parameters:
            self.parameters['magmoms'] = [atom.magmom for atom in self.atoms_sorted]

        # DFT + U uses a dictionary of: {symbol: {'L': int, 'U':float 'J':float}}
        # It is incompatible with setups because the dictionary doesn't allow
        # differentiating between symbols
        if 'ldau_luj' in self.parameters:
            if 'setups' in self.parameters:
                raise Exception('setups and ldau_luj is not supported.')

            atom_types = [x[0] if isinstance(x[0], str)
                          else self.atoms[x[0]].symbol
                          for x in self.ppp_list]

            d = self.parameters['ldau_luj']
            self.parameters['ldaul'] = [d[sym]['L'] for sym in atom_types]
            self.parameters['ldauu'] = [d[sym]['U'] for sym in atom_types]
            self.parameters['ldauj'] = [d[sym]['J'] for sym in atom_types]

    def set_label(self, label):
        """Set working directory.

        In VASP there is no prefix, only the working directory.

        """

        if label is None:
            self.directory = "."
            self.prefix = None
        else:
            self.directory, self.prefix = label, None

    def check_state(self, atoms):
        system_changes = FileIOCalculator.check_state(self, atoms)
        # Ignore boundary conditions:
        if 'pbc' in system_changes:
            system_changes.remove('pbc')

        return system_changes

    def set(self, **kwargs):
        """Set parameters with keyword=value pairs.

        calc.set(xc='PBE')

        TODO: I believe this will be problematic for special keywords

        """
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def read_results(self):
        """Read energy, forces, ... from output file(s)."""
        from ase.io.vasp import read_vasp_xml

        atoms = read_vasp_xml(os.path.join(self.directory,
                                           'vasprun.xml')).next()

        self.results['energy'] = atoms.get_potential_energy()
        self.results['forces'] = atoms.get_forces()

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=None):
        """Runs a calculation.

        No checking is currently done to see if one is necessary.

        """

        # This must setup some things.
        Calculator.calculate(self, atoms, properties, system_changes)

        self.write_input(self.atoms, properties, system_changes)
        if self.command is None:
            raise RuntimeError('Please set $%s environment variable ' %
                               ('ASE_' + self.name.upper() + '_COMMAND') +
                               'or supply the command keyword')

        olddir = os.getcwd()
        try:
            os.chdir(self.directory)
            errorcode = subprocess.call(self.command, shell=True)
        finally:
            os.chdir(olddir)

        if errorcode:
            raise RuntimeError('%s returned an error: %d' %
                               (self.name, errorcode))
        self.read_results()
