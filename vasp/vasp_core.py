"""The core Vasp calculator.

I aim to keep this file at a minimum. Hence, many logically grouped
class methods are actually imported at the end.

"""

import os
import subprocess
import numpy as np
from ase.calculators.calculator import Calculator
from ase.calculators.calculator import FileIOCalculator

import exceptions
from vasp import log


def VaspExceptionHandler(exc_type, exc_value, exc_traceback):
    """Handle exceptions."""
    if exc_type == exceptions.VaspSubmitted:
        print exc_value
        return None
    elif exc_type == exceptions.VaspQueued:
        print exc_value
        return None
    elif exc_type == KeyError:
        # This is a common error getting things from a dictionary,
        # especially the results dictionary.
        return None
    raise


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
    debug = None

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
        self.debug = debug
        # Here we decorate all the methods with the tryit decorator
        # unless we are debugging.
        if debug is not None:
            log.setLevel(debug)
        else:
            from vasp import tryit
            for attr in self.__dict__:
                if callable(getattr(self, attr)):
                    setattr(self, attr, tryit(getattr(self, attr)))

            for attr in Calculator.__dict__:
                if callable(getattr(Calculator, attr)):
                    setattr(Calculator, attr, tryit(getattr(Calculator, attr)))

            for attr in FileIOCalculator.__dict__:
                if callable(getattr(FileIOCalculator, attr)):
                    setattr(FileIOCalculator, attr,
                            tryit(getattr(FileIOCalculator, attr)))

        # This makes sure directories are created. We do not pass
        # kwargs here. Some of the special kwargs cannot be set at
        # this point.
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms)

        self.exception_handler = exception_handler

        # Ispin is especially problematic. In the set_ispin_dict
        # function I automatically set magmoms, but this requires
        # access to a sorted atoms, which we do not have yet. We save
        # it here and deal with it later.

        if 'ispin' in kwargs:
            ispin = kwargs['ispin']
            del kwargs['ispin']
        else:
            ispin = None

        if 'ldau_luj' in kwargs:
            ldau_luj = kwargs['ldau_luj']
            del kwargs['ldau_luj']
        else:
            ldau_luj = None

        self.set(**kwargs)

        self.read_metadata()  # this is ok if METADATA does not exist.
        log.debug('init metadata = {}'.format(self.metadata))

        if atoms is not None:
            # Resort atoms and generate list of POTCARS/counts we
            # need.  This code relies on self.parameters because of pp
            # and setups.
            self.sort_atoms(atoms)
        else:
            import ase.io
            contcar = os.path.join(self.directory, 'CONTCAR')
            empty_contcar = False
            if os.path.exists(contcar):
                # make sure the contcar is not empty
                with open(contcar) as f:
                    if f.read() == '':
                        empty_contcar = True

            poscar = os.path.join(self.directory, 'POSCAR')

            if os.path.exists(contcar) and not empty_contcar:
                atoms = ase.io.read(contcar)[self.metadata['resort']]
            elif os.path.exists(poscar):
                atoms = ase.io.read(poscar)[self.metadata['resort']]
            self.sort_atoms(atoms)

        # We need to handle some special kwargs here. The idea is that
        # self.parameters is ready to be written out later. This has
        # to be done after the sorting is known so that arrays of
        # properties can be made the right size and in the right order

        # Spin-polarization.
        if ispin:
            d = self.set_ispin_dict(ispin)
            self.parameters.update(d)
        if ldau_luj:
            d = self.set_ldau_luj_dict(ldau_luj)
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
            d = self.xc_defaults[self.parameters['xc'].lower()]
            self.parameters.update(d)

        # Done with initialization from kwargs. Now we need to figure
        # out what do about existing calculations that may exist.
        state = self.get_state()

        if state in [Vasp.NEW, Vasp.EMPTY]:
            log.debug('Calculation is empty or new')
            pass

        elif state == Vasp.QUEUED:
            log.debug('Calculation is queued')
            atoms, params = self.read()
            if self.atoms is None:
                self.atoms = atoms
            self.parameters.update(params)
            # no need to read results here. They don't exist yet.

        elif state == Vasp.FINISHED:
            log.debug('Finished calculation.')
            atoms, params = self.read()
            self.parameters.update(params)
            self.read_results()

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
        return atoms[self.resort]

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

    def check_state(self, atoms=None):
        """Check if any changes exist that require new calculations."""
        if atoms is None:
            atoms = self.get_atoms()

        system_changes = FileIOCalculator.check_state(self, atoms)
        # Ignore boundary conditions:
        if 'pbc' in system_changes:
            system_changes.remove('pbc')

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
                file_params['xc'] = ex
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

        if not self.parameters == file_params:
            new_keys = set(self.parameters.keys()) - set(file_params.keys())
            missing_keys = (set(file_params.keys()) -
                            set(self.parameters.keys()))

            log.debug('New keys: {}'.format(new_keys))
            log.debug('Missing keys: {}'.format(missing_keys))
            system_changes += ['params_on_file']

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

        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()  # needs to be resorted
        stress = atoms.get_stress()

        if self.atoms is None:
            atoms = atoms[self.metadata['resort']]
            self.sort_atoms(atoms)
            self.atoms.set_calculator(self)
        else:
            # update the atoms
            self.atoms.positions = atoms.positions[self.metadata['resort']]
            self.atoms.cell = atoms.cell

        self.results['energy'] = energy
        self.results['forces'] = forces[self.resort]
        self.results['stress'] = stress

        magnetic_moment = 0
        magnetic_moments = np.zeros(len(atoms))
        if self.parameters.get('ispin', 0) == 2:
            lines = open(os.path.join(self.directory, 'OUTCAR'), 'r').readlines()
            for n, line in enumerate(lines):
                if line.rfind('number of electron  ') > -1:
                    magnetic_moment = float(line.split()[-1])

                if line.rfind('magnetization (x)') > -1:
                    for m in range(len(atoms)):
                        magnetic_moments[m] = float(lines[n + m + 4].split()[4])

            self.results['magmom'] = magnetic_moment
            self.results['magmoms'] = np.array(magnetic_moments)[self.resort]

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

    def abort(self):
        """Abort and exit the program the calculator is running in."""
        import sys
        sys.exit()

    def stop_if(self, condition=None):
        """Stop program if condition is truthy."""
        if condition:
            import sys
            sys.exit()

    def clone(self, newdir):
        """Copy the calculation directory to newdir and set label to
        newdir.

        """
        state = self.get_state()

        import shutil
        if not os.path.isdir(newdir):
            shutil.copytree(self.directory, newdir)

            # need some cleanup here. do not copy jobids, etc...
            md = os.path.join(newdir, 'METADATA')
            if os.path.exists(md):
                import json
                with open(md) as f:
                    j = json.loads(f.read())
                    if 'jobid' in j:
                        del j['jobid']
                with open(md, 'wb') as f:
                    f.write(json.dumps(j))

            # What survives depends on the state
            # delete these files if not finished.
            if state in [Vasp.QUEUED, Vasp.NOTFINISHED]:
                os.unlink(os.path.join(newdir, 'OUTCAR'))
                os.unlink(os.path.join(newdir, 'vasprun.xml'))

            if state in [Vasp.EMPTYCONTCAR]:
                os.unlink(os.path.join(newdir, 'OUTCAR'))
                os.unlink(os.path.join(newdir, 'vasprun.xml'))
                os.unlink(os.path.join(newdir, 'CONTCAR'))

        self.__init__(newdir)

    def get_state(self):
        """Determine calculation state based on directory contents.

        Returns an integer for the state.

        """

        base_input = [os.path.exists(os.path.join(self.directory, f))
                     for f in ['INCAR', 'POSCAR', 'POTCAR', 'KPOINTS']]
        
        # Some input does not exist
        if False in base_input:
            # some input file is missing
            return Vasp.EMPTY

        # Input files exist, but no jobid, and no output
        if (np.array(base_input).all()
            and 'jobid' not in self.metadata
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

        # Not in queue, and not finished, with empty contcar
        if not self.in_queue():
            if os.path.exists(self.contcar):
                with open(self.contcar) as f:
                    if f.read() == '':
                        return Vasp.EMPTYCONTCAR
        
        return Vasp.UNKNOWN
