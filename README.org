#+TITLE: A new ASE interface to Vasp
#+date: May 4, 2016

* Why?
To make it compliant with the FileIOCalculator in ase, and hopefully simplify it.

* Goals
1. Provide an interface that allows all vasp INCAR tags without
   further modification, i.e. they are mostly just kwargs to the
   calculator.
2. Eliminate need for long lists of supported keywords.
2. Provide  an interface that can be extended to allow validation, and
   user-defined special keywords.
3. Better treatment of the xc keyword to make it easier to use
   different functionals.
4. Smart restarts. We don't want to run unnecessary calculations if
   they are already done. We want one script that runs calculations
   and does analysis.

* Installation

Clone this repo somewhere. Add it to your PYTHONPATH and PATH. Set VASP_PP_PATH to point to the directory containing your POTCARs.

#+BEGIN_SRC python
export PYTHONPATH=~/kitchin-python/vasp:$PYTHONPATH
export PATH=~/kitchin-python/vasp/bin:$PATH

export VASP_PP_PATH=/opt/kitchingroup/vasp-5.3.5
#+END_SRC

These directories should contain your POTCARs.
#+BEGIN_SRC sh
ls -d $VASP_PP_PATH/potpaw_???
#+END_SRC

#+RESULTS:
: /opt/kitchingroup/vasp-5.3.5/potpaw_GGA
: /opt/kitchingroup/vasp-5.3.5/potpaw_LDA
: /opt/kitchingroup/vasp-5.3.5/potpaw_PBE

* Differences from ase.calculators.vasp
Most things are the same. Here are few differences.

** label is a directory and the first argument
"label" is the first argument to the calculator, and it specifies the directory where the results are. Almost all file-io is done by path, so few directory changes ever occur.

** perpetual restart mode
Always starts in "restart" mode. On initialization the calculator always updates from the file system first, including updating the atoms, then from arguments. This allows you to write one script to setup, run and perform analysis.

** special setup syntax
Special setups are now specified by a list of [atom_symbol, potcar_suffix]

In this example we use the potpaw_PBE/O_s/POTCAR.

#+BEGIN_SRC python
calc = Vasp('molecules/O_s',
          encut=300,
            xc='PBE',
            ispin=2,
            ismear=0,
            sigma=0.001,
            setups=[['O', '_s']], # specifies O_s potential
            atoms=atoms)
#+END_SRC

This was changed to help make resorting simpler and reliable.

** new rwigs syntax
 rwigs is now a dictionary of {atom-symbol: radius}. This makes it easier to correctly generate the INCAR.

** ADOS is part of Vasp
The syntax to get the 's' orbital on the 0-indexed atom is:
#+BEGIN_SRC python
energies, c_s = calc.get_ados(0, 's')
#+END_SRC

Only 's', 'p', and 'd' are currently supported.
** Integrated visualization
This will show you the trajectory of the geometry relaxation.
#+BEGIN_SRC python
from vasp import Vasp
calc = Vasp('molecules/h2o-relax-centered')
calc.view()
#+END_SRC

** New default parameters
These may change. We don't usually write out the charge and wavecar files because they are large. An exception is if nsw>0, then we do write out the wavecar file to facilitate restarts.

#+BEGIN_SRC python
from vasp import Vasp
print(Vasp.default_parameters)
#+END_SRC

#+RESULTS:
: {'lcharg': False, 'kpts': [1, 1, 1], 'ismear': 1, 'xc': 'PBE', 'lwave': False, 'sigma': 0.1, 'pp': 'PBE'}

** Automatic job submission and job management.
Calculations are automatically submitted to a queue system with well-defined exceptions to provide job management. The setup is somewhat general, and must be tuned for specific clusters.

** Built-in exception handling.
All functions are wrapped in exception handling code to make some things easy to handle.

** "Smart" kwarg expansion.
Some kwargs are special, e.g. you can set ispin=2 and the calculator will automatically set the magmom key from the atoms object.

** Native support for the ase-db.
We actually use the ase-db to store calculation information.

#+BEGIN_SRC python
from vasp import Vasp
from ase.db import connect

bond_lengths = [1.05, 1.1, 1.15, 1.2, 1.25]
calcs = [Vasp('molecules/co-{0}'.format(d)) for d in bond_lengths]

con = connect('co-database.db', append=False)
for atoms in [calc.get_atoms() for calc in calcs]:
    con.write(atoms)
#+END_SRC

** Validation of some kwargs.
The vasp.validate file defines validation functions for many keywords, as well as brief documentation for them. This is integrated with Emacs to provide tooltips and easy access to documentation while working.

** VASPRC
This is a configuration file that allows customization of how jobs are submitted and whether validation is performed.

* Examples of usage
** A simple CO calculation
This is the prototypical simple calculation.

#+BEGIN_SRC python
from ase import Atoms, Atom
from vasp import Vasp
from vasp.vasprc import VASPRC
VASPRC['mode'] = 'run'

co = Atoms([Atom('C', [0, 0, 0]),
            Atom('O', [1.2, 0, 0])],
           cell=(6., 6., 6.))

calc = Vasp('~/dft-book-new-vasp/molecules/simple-co',  # output dir
            xc='pbe',    # the exchange-correlation functional
            nbands=6,    # number of bands
            encut=350,   # planewave cutoff
            ismear=1,    # Methfessel-Paxton smearing
            sigma=0.01,  # very small smearing factor for a molecule
            atoms=co)

print('energy = {0} eV'.format(co.get_potential_energy()))
print(co.get_forces())
#+END_SRC

#+RESULTS:
: energy = -14.69111507 eV
: [[ 5.09138064  0.          0.        ]
:  [-5.09138064  0.          0.        ]]

** A functional approach to calculations

Here we use list comprehensions to calculate the energy as a function of bond lengths.
#+BEGIN_SRC python :results output :exports both
from vasp import Vasp
from ase import Atom, Atoms
import logging

bond_lengths = [1.05, 1.1, 1.15, 1.2, 1.25]

ATOMS = [Atoms([Atom('C', [0, 0, 0]),
                Atom('O', [d, 0, 0])],
               cell=(6, 6, 6))
         for d in bond_lengths]

calcs = [Vasp('~/dft-book-new-vasp/molecules/co-{0}'.format(d),  # output dir
                xc='PBE',
                nbands=6,
                encut=350,
                ismear=1,
                sigma=0.01, debug=True,
                atoms=atoms)
         for d, atoms in zip(bond_lengths, ATOMS)]

energies = [atoms.get_potential_energy() for atoms in ATOMS]

print(energies)
#+END_SRC
 tpptree
#+RESULTS:
: [-14.21584765, -14.72174343, -14.84115208, -14.69111507, -14.35508371]

** Some new ideas in job management
By default, many exceptions are handled automatically, and if calculations are not finished the quantities are returned as None. This leads to some challenges if you want to do analysis before the results are ready.

Our workflow relies on asynchronously running jobs in a queue. To avoid blocking scripts, we setup everything so that scripts just exit if they cannot continue, and we rerun them later.

We provide the following tools for handling these situations.

*** calc.abort()
The abort function simply exits the program when called.
#+BEGIN_SRC python
from vasp import Vasp
from ase.lattice.cubic import BodyCenteredCubic

atoms = BodyCenteredCubic(directions=[[1, 0, 0],
                                      [0, 1, 0],
                                      [0, 0, 1]],
                                      size=(1, 1, 1),
                                      symbol='Fe')

NUPDOWNS = [0.0, 2.0, 4.0, 5.0, 6.0, 8.0]
energies = []
for B in NUPDOWNS:
    calc = Vasp('bulk/Fe-bcc-fixedmagmom-{0:1.2f}'.format(B),
                xc='PBE',
                encut=300,
                kpts=(4, 4, 4),
                ispin=2,
                nupdown=B,
                atoms=atoms)
    energies.append(atoms.get_potential_energy())

if None in energies:
    calc.abort()

# some analysis that depends on all energies being present
#+END_SRC

*** calc.wait()
The wait function does not actually wait. It does try to get the energy and run the job, and if it is not ready, it exits. The name or action of this function may change.

#+BEGIN_SRC python
from vasp import Vasp
from ase.lattice.cubic import FaceCenteredCubic

atoms = FaceCenteredCubic(symbol='Al')

calc = Vasp('bulk/Al-bulk',
            xc='PBE',
            kpts=(12, 12, 12),
            encut=350,
            prec='High',
            isif=3,
            nsw=30,
            ibrion=1,
            atoms=atoms)
calc.wait()

# some analysis that depends on the calculation being done
#+END_SRC

*** calc.stop_if(condition)
Sometimes you would like some condition to determine if you stop. This is a one line version of the if statement combined with calc.abort()

#+BEGIN_SRC python
from vasp import Vasp
from ase import Atom, Atoms
import numpy as np
# fcc
LC = [3.5, 3.55, 3.6, 3.65, 3.7, 3.75]
volumes, energies = [], []
for a in LC:
    atoms = Atoms([Atom('Ni', (0, 0, 0), magmom=2.5)],
                  cell=0.5 * a * np.array([[1.0, 1.0, 0.0],
                                           [0.0, 1.0, 1.0],
                                           [1.0, 0.0, 1.0]]))

    calc = Vasp('bulk/Ni-{0}'.format(a),
                xc='PBE',
                encut=350,
                kpts=(12, 12, 12),
                ispin=2,
                atoms=atoms)
    energies.append(calc.potential_energy)
    volumes.append(atoms.get_volume())

calc.stop_if(None in energies)

# some analysis requireing all the energies.
#+END_SRC

** Run simulations with a Lisp
One of my motivations for the rewrite was to enable me to use Hy (http://docs.hylang.org/en/latest/) in these simulations. Hy is a Lisp that runs Python. Here is an example calculation. This might be interesting because it allows you to write macros. I am not sure what I will do that yet, but I look forward to trying it out.

#+BEGIN_SRC hy
(import [ase [Atom Atoms]])
(import [vasp [Vasp]])

(setv co (Atoms [(Atom "C" [0.0 0.0 0.0])
                 (Atom "O" [1.2 0.0 0.0])]
                :cell [6.0 6.0 6.0]))

(setv calc (Vasp "~/dft-book-new-vasp/molecules/simple-co-hy"
                 :xc "pbe"
                 :nbands 6
                 :encut 350
                 :ismear 1
                 :sigma 0.01
                 :atoms co))

(print (.format "energy = {0} eV"
                (.get_potential_energy co)))

(print calc.potential_energy)
(print (.get_forces co))
#+END_SRC

#+RESULTS:
: energy = -14.69111507 eV
: -14.69111507
: [[ 5.09138064  0.          0.        ]
:  [-5.09138064  0.          0.        ]]

** vaspsum
This command line utility provides a variety of ways to summarize a calculation. For example, you can use this to print the input files:

#+BEGIN_SRC sh
vaspsum --vasp ~/dft-book-new-vasp/molecules/simple-co
#+END_SRC

#+RESULTS:
#+begin_example
INCAR
-----
INCAR created by Atomic Simulation Environment
 ENCUT = 350
 LCHARG = .FALSE.
 NBANDS = 6
 ISMEAR = 1
 LWAVE = .FALSE.
 SIGMA = 0.01


POSCAR
------
 C  O
 1.0000000000000000
     6.0000000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000    6.0000000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000    6.0000000000000000
   1   1
Cartesian
  0.0000000000000000  0.0000000000000000  0.0000000000000000
  1.2000000000000000  0.0000000000000000  0.0000000000000000


KPOINTS
-------
KPOINTS created by Atomic Simulation Environment
0
Monkhorst-Pack
1 1 1
0.0 0.0 0.0


POTCAR
------
cat $VASP_PP_PATH/potpaw_PBE/C/POTCAR $VASP_PP_PATH/potpaw_PBE/O/POTCAR > POTCAR
#+end_example

Or this to output the ase-db json.
#+BEGIN_SRC sh
vaspsum --json ~/dft-book-new-vasp/molecules/simple-co
#+END_SRC

#+RESULTS:
#+begin_example
json:  {'lcharg': False, 'pp': 'PBE', 'nbands': 6, 'xc': 'pbe', 'ismear': 1, 'lwave': False, 'sigma': 0.01, 'kpts': [1, 1, 1], 'encut': 350}
{"1": {
 "calculator": "vasp",
 "calculator_parameters": {"xc": "pbe", "nbands": 6, "sigma": 0.01, "encut": 350},
 "cell": [[6.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 6.0]],
 "charges": [null, null],
 "ctime": 16.380341757550546,
 "data": {"resort": [0, 1], "ppp_list": [["C", "potpaw_PBE/C/POTCAR", 1], ["O", "potpaw_PBE/O/POTCAR", 1]], "parameters": {"lcharg": false, "pp": "PBE", "nbands": 6, "xc": "pbe", "ismear": 1, "lwave": false, "sigma": 0.01, "kpts": [1, 1, 1], "encut": 350}},
 "energy": -14.69111507,
 "forces": [[5.09138064, 0.0, 0.0], [-5.09138064, 0.0, 0.0]],
 "key_value_pairs": {"path": "/home-research/jkitchin/dft-book-new-vasp/molecules/simple-co"},
 "magmom": 0,
 "magmoms": [0.0, 0.0],
 "mtime": 16.380341757550546,
 "numbers": [6, 8],
 "pbc": [true, true, true],
 "positions": [[0.0, 0.0, 0.0], [1.2000000000000002, 0.0, 0.0]],
 "stress": [0.041455596684986905, 0.01094970637584278, 0.01094970637584278, -0.0, -0.0, -0.0],
 "unique_id": "671032550621923e208be983ce744d24",
 "user": "jkitchin"},
"ids": [1],
"nextid": 2}

#+end_example

* vaspy-mode
We provide vaspy-mode to enhance using vasp in Emacs. The main feature it provides is syntax highlighting on vasp keywords with a tooltip on them showing the first line of the validation docstring, and making them clickable to show the whole docstring.

Add this to your Emacs initialization file (obviously change the path to where you installed the vasp module.

#+BEGIN_SRC emacs-lisp
(add-to-list 'load-path "~/kitchin-python/vasp")
(require 'vaspy-mode)
#+END_SRC

#+RESULTS:
: vaspy-mode
* Documentation
Here is a list of commands, their docstrings, links to the code and the code for reference.
#+BEGIN_SRC python :results raw
import inspect

from vasp.vasprc import VASPRC
VASPRC['handle_exceptions'] = False

from vasp import Vasp

print "** Vasp functions"

for attr in sorted(Vasp.__dict__.keys()):
    if callable(getattr(Vasp, attr)):
        argspec = inspect.getargspec(Vasp.__dict__[attr])
        args = argspec.args
        varargs = argspec.varargs
        kwargs = argspec.keywords


        defaults = argspec.defaults
        if defaults is not None:
            N = len(args) - len(defaults)

            argstring = ', '.join(args[0: N])
            argstring += ', ' + ', '.join(['{}={}'.format(a, b)
                                           for a, b in zip(args[N:], defaults)])
        else:
            argstring = ', '.join(args)

        if varargs is not None:
            argstring += ', *{}'.format(varargs)

        if kwargs is not None:
            argstring += ', **{}'.format(kwargs)
        pyfile = inspect.getfile(getattr(Vasp, attr))
        source, lineno = inspect.getsourcelines(getattr(Vasp, attr))
        print('*** Vasp.{0}\nargs = ({1})\n\n{2}\n'.format(attr,
                                                      argstring,
                                                      Vasp.__dict__[attr].__doc__))

        print '[[./{}::{}]]'.format(pyfile, lineno)

        print """#+BEGIN_SRC python
{}
,#+END_SRC
""".format(''.join(source))
#+END_SRC

#+RESULTS:
** Vasp functions
*** Vasp.__init__
args = (self, label, restart=True, ignore_bad_restart_file=False, atoms=None, scratch=None, debug=None, exception_handler=<function VaspExceptionHandler at 0x1a88ae60>, **kwargs)

Create a Vasp calculator.

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



[[./vasp/vasp_core.py::137]]
#+BEGIN_SRC python
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
        # cast as str in case label is unicode, i.e. if it is from hy.
        self.set_label(label)
        self.debug = debug
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
                a.pbs = [True, True, True]
            self.neb = True

        # We do not pass kwargs here. Some of the special kwargs
        # cannot be set at this point since they need to know about
        # the atoms and parameters. This reads params and results from
        # existing files if they are there. It calls self.read(). It
        # should update the atoms from what is on file.

        if self.neb is not None:
            FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                      str(label))
            self.neb = atoms
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

        # These depend on having atoms already.
        if ispin is not None:
            self.set(**self.set_ispin_dict(ispin))

        if rwigs is not None:
            self.set(**self.set_rwigs_dict(rwigs))

        # Finally run validate functions
        if VASPRC['validate']:
            for key, val in self.parameters.iteritems():
                if key in validate.__dict__:
                    f = validate.__dict__[key]
                    f(self, val)
                else:
                    warnings.warn('No validation for {}'.format(key))

#+END_SRC

*** Vasp.__str__
args = (self)

Pretty representation of a calculation.

        TODO: make more like jaspsum.



[[./vasp/vasp_core.py::361]]
#+BEGIN_SRC python
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

#+END_SRC

*** Vasp.abort
args = (self)

Abort and exit the program the calculator is running in.

[[./vasp/vasp_core.py::548]]
#+BEGIN_SRC python
    def abort(self):
        """Abort and exit the program the calculator is running in."""
        import sys
        sys.exit()

#+END_SRC

*** Vasp.attach_charges
args = (self, fileobj=None, displacement=0.0001)

Attach the charges from the fileobj to the atoms on the calculator.

    This is a modified version of the attach_charges function in
    ase.io.bader to work better with VASP.
    Does not require the atom positions to be in Bohr and references
    the charge to the ZVAL in the POTCAR


Monkey-patch defined in vasp/bader.py at line 9

[[./vasp/bader.py::9]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def attach_charges(self, fileobj=None, displacement=1e-4):
    """Attach the charges from the fileobj to the atoms on the calculator.

    This is a modified version of the attach_charges function in
    ase.io.bader to work better with VASP.
    Does not require the atom positions to be in Bohr and references
    the charge to the ZVAL in the POTCAR
    """
    if fileobj is None:
        fileobj = os.path.join(self.directory, 'ACF.dat')

    if isinstance(fileobj, str):
        fileobj = open(fileobj)
        f_open = True

    # First get a dictionary of ZVALS from the pseudopotentials
    LOP = self.get_pseudopotentials()
    ppp = os.environ['VASP_PP_PATH']

    zval = {}
    for sym, ppath, hash in LOP:
        fullpath = os.path.join(ppp, ppath)
        z = get_ZVAL(fullpath)
        zval[sym] = z

    atoms = self.atoms
    # Get sorted symbols and positions according to POSCAR and ACF.dat
    symbols = np.array(atoms.get_chemical_symbols())[self.resort]
    positions = atoms.get_positions()[self.resort]

    charges = []
    sep = '---------------'
    i = 0  # Counter for the lines
    k = 0  # Counter of sep
    assume6columns = False
    for line in fileobj:
        if line[0] == '\n':  # check if there is an empty line in the
            i -= 1           # head of ACF.dat file
        if i == 0:
            headings = line
            if 'BADER' in headings.split():
                j = headings.split().index('BADER')
            elif 'CHARGE' in headings.split():
                j = headings.split().index('CHARGE')
            else:
                print('Can\'t find keyword "BADER" or "CHARGE".'
                      ' Assuming the ACF.dat file has 6 columns.')
                j = 4
                assume6columns = True
        if sep in line:  # Stop at last seperator line
            if k == 1:
                break
            k += 1
        if not i > 1:
            pass
        else:
            words = line.split()
            if assume6columns is True:
                if len(words) != 6:
                    raise IOError('Number of columns in ACF file incorrect!\n'
                                  'Check that Bader program version >= 0.25')

            sym = symbols[int(words[0]) - 1]
            charges.append(zval[sym] - float(words[j]))

            if displacement is not None:
                # check if the atom positions match
                xyz = np.array([float(w) for w in words[1:4]])
                assert (np.linalg.norm(positions[int(words[0]) - 1] - xyz)
                        < displacement)
        i += 1

    if f_open:
        fileobj.close()

    # Now attach the resorted charges to the atom
    charges = np.array(charges)[self.resort]
    for atom in self.atoms:
        atom.charge = charges[atom.index]

#+END_SRC

*** Vasp.bader
args = (self, cmd=None, ref=False, verbose=False, overwrite=False)

Performs bader analysis for a calculation.
    Follows defaults unless full shell command is specified
    Does not overwrite existing files if overwrite=False
    If ref = True, tries to reference the charge density to
    the sum of AECCAR0 and AECCAR2
    Requires the bader.pl (and chgsum.pl) script to be in the system PATH


Monkey-patch defined in vasp/bader.py at line 108

[[./vasp/bader.py::108]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def bader(self, cmd=None, ref=False, verbose=False, overwrite=False):
    """Performs bader analysis for a calculation.
    Follows defaults unless full shell command is specified
    Does not overwrite existing files if overwrite=False
    If ref = True, tries to reference the charge density to
    the sum of AECCAR0 and AECCAR2
    Requires the bader.pl (and chgsum.pl) script to be in the system PATH
    """
    cwd = os.getcwd()
    try:
        os.chdir(self.directory)

        if 'ACF.dat' in os.listdir(".") and not overwrite:
            self.attach_charges()
            return

        if cmd is None:
            if ref:
                self.chgsum()
                cmdlist = ['bader',
                           'CHGCAR',
                           '-ref',
                           'CHGCAR_sum']
            else:
                cmdlist = ['bader', 'CHGCAR']
        elif type(cmd) is str:
            cmdlist = cmd.split()
        elif type(cmd) is list:
            cmdlist = cmd

        p = Popen(cmdlist, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        if out == '' or err != '':
            raise Exception('Cannot perform Bader:\n\n{0}'.format(err))
        elif verbose:
            print('Bader completed for {0}'.format(self.vaspdir))

        self.attach_charges('ACF.dat')
    finally:
        os.chdir(cwd)

#+END_SRC

*** Vasp.calculate
args = (self, atoms=None, properties=['energy'], system_changes=None)

Monkey patch to submit job through the queue.
    If this is called, then the calculator thinks a job should be run.
    If we are in the queue, we should run it, otherwise, a job should
    be submitted.


Monkey-patch defined in vasp/runner.py at line 61

[[./vasp/runner.py::61]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def calculate(self, atoms=None, properties=['energy'],
              system_changes=None):
    """Monkey patch to submit job through the queue.
    If this is called, then the calculator thinks a job should be run.
    If we are in the queue, we should run it, otherwise, a job should
    be submitted.
    """
    log.debug('In queue: {}'.format(self.in_queue()))
    if self.in_queue():
        raise VaspQueued('{} Queued: {}'.format(self.directory,
                                                self.get_db('jobid')))

    # not in queue. Delete the jobid
    if self.get_db('jobid') is not None:
        self.write_db(jobid=None)

        # we should check for errors here.
        self.read_results()
        return

    if (not self.calculation_required(atoms, ['energy'])
        and not self.check_state()):
        print('No calculation_required.')
        self.read_results()
        return

    # The subclass implementation should first call this
    # implementation to set the atoms attribute.
    Calculator.calculate(self, atoms, properties, system_changes)

    self.write_input(atoms, properties, system_changes)
    if self.parameters.get('luse_vdw', False):
        kernel = os.path.join(self.directory, 'vdw_kernel.bindat')
        if not os.path.exists(kernel):
            os.symlink(VASPRC['vdw_kernel.bindat'], kernel)

    # if we are in the queue and vasp is called or if we want to use
    # mode='run' , we should just run the job. First, we consider how.
    if 'PBS_O_WORKDIR' in os.environ or VASPRC['mode'] == 'run':
        if 'PBS_NODEFILE' in os.environ:
            # we are in the queue. determine if we should run serial
            # or parallel
            NPROCS = len(open(os.environ['PBS_NODEFILE']).readlines())
            log.debug('Found {0} PROCS'.format(NPROCS))
            if NPROCS == 1:
                # no question. running in serial.
                vaspcmd = VASPRC['vasp.executable.serial']
                log.debug('NPROCS = 1. running in serial')
                exitcode = os.system(vaspcmd)
                return exitcode
            else:
                # vanilla MPI run. multiprocessing does not work on more
                # than one node, and you must specify in VASPRC to use it
                if (VASPRC['queue.nodes'] > 1
                    or (VASPRC['queue.nodes'] == 1
                        and VASPRC['queue.ppn'] > 1
                        and (VASPRC['multiprocessing.cores_per_process']
                             == 'None'))):
                    s = 'queue.nodes = {0}'.format(VASPRC['queue.nodes'])
                    log.debug(s)
                    log.debug('queue.ppn = {0}'.format(VASPRC['queue.ppn']))
                    mpc = VASPRC['multiprocessing.cores_per_process']
                    log.debug('multiprocessing.cores_per_process'
                              '= {0}'.format(mpc))
                    log.debug('running vanilla MPI job')

                    log.debug('MPI NPROCS = {}'.format(NPROCS))
                    vaspcmd = VASPRC['vasp.executable.parallel']
                    parcmd = 'mpirun -np %i %s' % (NPROCS, vaspcmd)
                    exitcode = os.system(parcmd)
                    return exitcode
                else:
                    # we need to run an MPI job on cores_per_process
                    if VASPRC['multiprocessing.cores_per_process'] == 1:
                        log.debug('running single core multiprocessing job')
                        vaspcmd = VASPRC['vasp.executable.serial']
                        exitcode = os.system(vaspcmd)
                    elif VASPRC['multiprocessing.cores_per_process'] > 1:
                        log.debug('running mpi multiprocessing job')
                        NPROCS = VASPRC['multiprocessing.cores_per_process']

                        vaspcmd = VASPRC['vasp.executable.parallel']
                        parcmd = 'mpirun -np %i %s' % (NPROCS, vaspcmd)
                        exitcode = os.system(parcmd)
                        return exitcode
        else:
            # probably running at cmd line, in serial.
            try:
                cwd = os.getcwd()
                os.chdir(self.directory)
                vaspcmd = VASPRC['vasp.executable.serial']
                status, output, err = getstatusoutput(vaspcmd,
                                                      stdout=subprocess.PIPE,
                                                      stderr=subprocess.PIPE)
                if status == 0:
                    self.read_results()
                    return True
                else:
                    return output
            finally:
                os.chdir(cwd)
        # end

    # if you get here, a job is getting submitted
    CWD = os.getcwd()
    VASPDIR = self.directory
    script = """
#!/bin/bash
cd {CWD}  # this is the current working directory
cd {VASPDIR}  # this is the vasp directory
runvasp.py     # this is the vasp command
#end""".format(**locals())

    jobname = VASPDIR
    log.debug('{0} will be the jobname.'.format(jobname))
    log.debug('-l nodes={0}:ppn={1}'.format(VASPRC['queue.nodes'],
                                            VASPRC['queue.ppn']))

    cmdlist = ['{0}'.format(VASPRC['queue.command'])]
    cmdlist += ['-o', VASPDIR]
    cmdlist += [option for option in VASPRC['queue.options'].split()]
    cmdlist += ['-N', '{0}'.format(jobname),
                '-l walltime={0}'.format(VASPRC['queue.walltime']),
                '-l nodes={0}:ppn={1}'.format(VASPRC['queue.nodes'],
                                              VASPRC['queue.ppn']),
                '-l mem={0}'.format(VASPRC['queue.mem'])]
    log.debug('{0}'.format(' '.join(cmdlist)))
    p = subprocess.Popen(cmdlist,
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    log.debug(script)

    out, err = p.communicate(script)

    if out == '' or err != '':
        raise Exception('something went wrong in qsub:\n\n{0}'.format(err))

    self.write_db(jobid=out.strip())

    raise VaspSubmitted('{} submitted: {}'.format(self.directory,
                                                  out.strip()))

#+END_SRC

*** Vasp.calculation_required
args = (self, atoms=None, properties=['energy'])

Returns if a calculation is needed.

[[./vasp/vasp_core.py::491]]
#+BEGIN_SRC python
    def calculation_required(self, atoms=None, properties=['energy']):
        """Returns if a calculation is needed."""

        if atoms is None:
            atoms = self.get_atoms()

        system_changes = self.check_state(atoms)
        if system_changes:
            print('Calculation needed for {}'.format(system_changes))
            return True

        for name in properties:
            if name not in self.results:
                print('{} not in {}. Calc required.'.format(name,
                                                                self.results))
                return True

        # if the calculation is finished we do not need to run.
        if os.path.exists(self.outcar):
            with open(self.outcar) as f:
                lines = f.readlines()
                if 'Voluntary context switches:' in lines[-1]:
                    return False

#+END_SRC

*** Vasp.check_state
args = (self, atoms=None)

Check if any changes exist that require new calculations.

[[./vasp/vasp_core.py::406]]
#+BEGIN_SRC python
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

        if not self.parameters == file_params:
            new_keys = set(self.parameters.keys()) - set(file_params.keys())
            missing_keys = (set(file_params.keys()) -
                            set(self.parameters.keys()))
            log.debug('New keys: {}'.format(new_keys))
            log.debug('Missing keys: {}'.format(missing_keys))
            system_changes += ['params_on_file']

        return system_changes

#+END_SRC

*** Vasp.chgsum
args = (self)

Uses the chgsum.pl utility to sum over the AECCAR0 and AECCAR2 files.

Monkey-patch defined in vasp/bader.py at line 91

[[./vasp/bader.py::91]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def chgsum(self):
    """Uses the chgsum.pl utility to sum over the AECCAR0 and AECCAR2 files."""
    cwd = os.getcwd()
    try:
        os.chdir(self.directory)
        cmdlist = ['chgsum.pl',
                   'AECCAR0',
                   'AECCAR2']
        p = Popen(cmdlist, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        if out == '' or err != '':
            raise Exception('Cannot perform chgsum:\n\n{0}'.format(err))
    finally:
        os.chdir(cwd)

#+END_SRC

*** Vasp.clone
args = (self, newdir)

Copy the calculation directory to newdir and set label to
        newdir.



[[./vasp/vasp_core.py::567]]
#+BEGIN_SRC python
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

        self.__init__(newdir)
        self.write_db(jobid=None, path=newdir)

#+END_SRC

*** Vasp.describe
args = (self, long=False)

Describe the parameters used with docstrings in vasp.validate.

[[./vasp/vasp_core.py::721]]
#+BEGIN_SRC python
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

#+END_SRC

*** Vasp.get_ados
args = (self, atom_index, orbital, spin=1, efermi=None)

Return Atom projected DOS for atom index, orbital and spin.

    orbital: string ['s', 'p', 'd']

    If efermi is not None, use this value as 0.0.

    :returns: (energies, ados)



Monkey-patch defined in vasp/getters.py at line 182

[[./vasp/getters.py::182]]
#+BEGIN_SRC python
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

#+END_SRC

*** Vasp.get_beefens
args = (self, n=-1)

Get the BEEFens 2000 ensemble energies from the OUTCAR.
    This only works with Vasp 5.3.5 compiled with libbeef.
    I am pretty sure this array is the deviations from the total
    energy. There are usually 2000 of these, but it is not clear this will
    always be the case. I assume the 2000 entries are always in the same
    order, so you can calculate ensemble energy differences for reactions,
    as long as the number of samples in the ensemble is the same.
    There is usually more than one BEEFens section. By default we
    return the last one. Choose another one with the the :par: n.
    see http://suncat.slac.stanford.edu/facility/software/functional/


Monkey-patch defined in vasp/getters.py at line 39

[[./vasp/getters.py::39]]
#+BEGIN_SRC python
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

#+END_SRC

*** Vasp.get_charge_density
args = (self, spin=0, filename=None)

Returns x, y, and z coordinate and charge density arrays.

    Supported file formats: CHG, CHGCAR
    :param int spin: an integer
    :returns: x, y, z, charge density arrays
    :rtype: 3-d numpy arrays
    Relies on :func:`ase.calculators.vasp.VaspChargeDensity`.


Monkey-patch defined in vasp/getters.py at line 327

[[./vasp/getters.py::327]]
#+BEGIN_SRC python
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

#+END_SRC

*** Vasp.get_db
args = (self, *keys)

Retrieve values for each key in keys.

    First look for key/value, then in data.



Monkey-patch defined in vasp/getters.py at line 12

[[./vasp/getters.py::12]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def get_db(self, *keys):
    """Retrieve values for each key in keys.

    First look for key/value, then in data.

    """
    dbfile = os.path.join(self.directory, 'DB.db')

    if not os.path.exists(dbfile):
        return [None for key in keys] if len(keys) > 1 else None

    vals = [None for key in keys]
    from ase.db import connect

    with connect(dbfile) as con:
        try:
            at = con.get(id=1)
            for i, key in enumerate(keys):
                vals[i] = (at.key_value_pairs.get(key, None)
                           or at.data.get(key, None))
        except KeyError, e:
            if e.message == 'no match':
                pass
    return vals if len(vals) > 1 else vals[0]

#+END_SRC

*** Vasp.get_default_number_of_electrons
args = (self, filename=None)

Return the default electrons for each species.

Monkey-patch defined in vasp/getters.py at line 243

[[./vasp/getters.py::243]]
#+BEGIN_SRC python
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

#+END_SRC

*** Vasp.get_dipole_moment
args = (self, atoms=None)

Return dipole_moment.

    dipole_moment = ((dipole_vector**2).sum())**0.5/Debye



Monkey-patch defined in vasp/getters.py at line 468

[[./vasp/getters.py::468]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def get_dipole_moment(self, atoms=None):
    """Return dipole_moment.

    dipole_moment = ((dipole_vector**2).sum())**0.5/Debye

    """
    self.update()

    dv = self.get_dipole_vector(atoms)

    from ase.units import Debye
    return ((dv ** 2).sum()) ** 0.5 / Debye

#+END_SRC

*** Vasp.get_dipole_vector
args = (self, atoms=None)

Tries to return the dipole vector of the unit cell in atomic units.

    Returns None when CHG file is empty/not-present.



Monkey-patch defined in vasp/getters.py at line 405

[[./vasp/getters.py::405]]
#+BEGIN_SRC python
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
        fullpath = os.path.join(ppp, ppath)
        z = get_ZVAL(fullpath)
        zval[sym] = z

    ion_charge_center = np.array([0.0, 0.0, 0.0])
    total_ion_charge = 0.0
    for atom in atoms:
        Z = zval[atom.symbol]
        total_ion_charge += Z
        pos = atom.position
        ion_charge_center += Z * pos

    ion_charge_center /= total_ion_charge
    ion_dipole_moment = ion_charge_center * total_ion_charge

    dipole_vector = (ion_dipole_moment + electron_dipole_moment)

    return dipole_vector

#+END_SRC

*** Vasp.get_eigenvalues
args = (self, kpt=0, spin=1)

Return array of eigenvalues for kpt and spin.

Monkey-patch defined in vasp/getters.py at line 144

[[./vasp/getters.py::144]]
#+BEGIN_SRC python
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

#+END_SRC

*** Vasp.get_elapsed_time
args = (self)

Return elapsed calculation time in seconds from the OUTCAR file.

Monkey-patch defined in vasp/getters.py at line 223

[[./vasp/getters.py::223]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def get_elapsed_time(self):
    """Return elapsed calculation time in seconds from the OUTCAR file."""
    self.update()
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

#+END_SRC

*** Vasp.get_electron_density_center
args = (self, spin=0, scaled=True)

Returns center of electron density.
    If scaled, use scaled coordinates, otherwise use cartesian
    coordinates.


Monkey-patch defined in vasp/getters.py at line 377

[[./vasp/getters.py::377]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def get_electron_density_center(self, spin=0, scaled=True):
    """Returns center of electron density.
    If scaled, use scaled coordinates, otherwise use cartesian
    coordinates.
    """
    self.update()
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

#+END_SRC

*** Vasp.get_elf
args = (self)

Returns x, y, z and electron localization function arrays.

Monkey-patch defined in vasp/getters.py at line 364

[[./vasp/getters.py::364]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def get_elf(self):
    """Returns x, y, z and electron localization function arrays."""
    assert self.parameters.get('lelf', None) is True,\
        "lelf is not set to True!"

    self.update()
    fname = os.path.join(self.directory, 'ELFCAR')
    x, y, z, data = get_volumetric_data(self, filename=fname)
    atoms = self.get_atoms()
    return x, y, z, data[0] * atoms.get_volume()

#+END_SRC

*** Vasp.get_fermi_level
args = (self)

Return the Fermi level.

Monkey-patch defined in vasp/getters.py at line 167

[[./vasp/getters.py::167]]
#+BEGIN_SRC python
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

#+END_SRC

*** Vasp.get_ibz_k_points
args = (self)

Return the irreducible k-points.

Monkey-patch defined in vasp/getters.py at line 63

[[./vasp/getters.py::63]]
#+BEGIN_SRC python
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

#+END_SRC

*** Vasp.get_infrared_intensities
args = (self)

Calculate infrared intensities of vibrational modes.

    Returns an array of normalized intensities for each vibrational
    mode. You should have run the vibrational calculation already. This
    function does not run it for you.

    python translation of # A utility for calculating the vibrational
    intensities from VASP output (OUTCAR) # (C) David Karhanek,
    2011-03-25, ICIQ Tarragona, Spain (www.iciq.es)
    http://homepage.univie.ac.at/david.karhanek/downloads.html#Entry02


Monkey-patch defined in vasp/vib.py at line 194

[[./vasp/vib.py::194]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def get_infrared_intensities(self):
    """Calculate infrared intensities of vibrational modes.

    Returns an array of normalized intensities for each vibrational
    mode. You should have run the vibrational calculation already. This
    function does not run it for you.

    python translation of # A utility for calculating the vibrational
    intensities from VASP output (OUTCAR) # (C) David Karhanek,
    2011-03-25, ICIQ Tarragona, Spain (www.iciq.es)
    http://homepage.univie.ac.at/david.karhanek/downloads.html#Entry02
    """
    assert self.parameters.get('lepsilon', None) is True
    assert self.parameters.get('nwrite', 0) == 3
    assert self.parameters.get('ibrion', 0) == 7

    self.update()

    atoms = read(os.path.join(self.directory, 'POSCAR'), format='vasp')
    NIONS = len(atoms)
    BORN_NROWS = NIONS * 4 + 1

    with open(os.path.join(self.directory, 'OUTCAR'), 'r') as f:
        alltext = f.read()
        f.seek(0)
        alllines = f.readlines()
        f.close()

    if 'BORN' not in alltext:
        raise Exception('Born effective charges missing. '
                        'Did you use IBRION=7 or 8?')

    if 'Eigenvectors after division by SQRT(mass)' not in alltext:
        raise Exception('You must rerun with NWRITE=3 to get '
                        'sqrt(mass) weighted eigenvectors')

    # get the Born charges
    for i, line in enumerate(alllines):
        if 'BORN EFFECTIVE CHARGES' in line:
            break

    BORN_MATRICES = []
    i += 2  # skip a line
    for j in range(NIONS):
        BM = []
        i += 1  # skips the ion count line
        for k in range(3):
            line = alllines[i]
            fields = line.split()
            BM.append([float(x) for x in fields[1:4]])
            i += 1  # advance a line
        BORN_MATRICES.append(BM)

    BORN_MATRICES = np.array(BORN_MATRICES)

    # Get the eigenvectors and eigenvalues.  maybe I can replace this
    # code with my other code. for now I just reproduce the count
    # number of vibs. this gets the number from outcar. it seems like
    # it should be known in advance unless constraints make it hard to
    # tell.

    # the next code in the shell script just copies code to eigenvectors.txt
    for i, line in enumerate(alllines):
        if 'Eigenvectors after division by SQRT(mass)' in line:
            break

    EIG_NVIBS = 0
    for line in alllines[i:]:
        if ('f' in line
            and 'THz' in line
            and 'cm-1' in line):
            EIG_NVIBS += 1

    EIG_NIONS = BORN_NROWS
    # I guess this counts blank rows and non-data rows
    # EIG_NROWS = (EIG_NIONS + 3) * EIG_NVIBS + 3

    # i is where the data starts
    i += 6

    EIGENVALUES = []
    EIGENVECTORS = []
    for j in range(EIG_NVIBS):
        mode = []
        EIGENVALUES.append(alllines[i])  # frequencies are here

        i += 1  # skip the frequency line
        i += 1  # skip the xyz line
        for k in range(3):
            fields = [float(x) for x in alllines[i].split()]
            mode.append(fields[3:])
            i += 1
        EIGENVECTORS.append(mode)
        i += 1  # skip blank line

    EIGENVECTORS = np.array(EIGENVECTORS)

    # now we are ready to compute intensities. see
    # http://othes.univie.ac.at/10117/1/2010-05-05_0547640.pdf, page
    # 21.

    # I(\omega) = \sum_{\alpha=1}^3 |
    # \sum_{l=1}^M \sum_{\beta=1}^3 Z_{\alpha\beta}(l)e_{\beta}(l)|^2

    # omega is the vibrational mode
    # alpha, beta are the cartesian polarizations
    # l is the atom number
    # e_beta is the eigenvector of the mode

    intensities = []

    for mode in range(len(EIGENVECTORS)):
        S = 0  # This is the triple sum
        for alpha in [0, 1, 2]:
            s = 0
            for l in [0, 1, 2]:  # this is the atom number
                for beta in [0, 1, 2]:
                    e = EIGENVECTORS[mode][l]
                    Zab = BORN_MATRICES[l][alpha][beta]

                    s += Zab * e[beta]
            S += s ** 2
        intensities.append(S)

    intensities = np.array(intensities) / max(intensities)
    return intensities

#+END_SRC

*** Vasp.get_k_point_weights
args = (self)

Return the k-point weights.

Monkey-patch defined in vasp/getters.py at line 118

[[./vasp/getters.py::118]]
#+BEGIN_SRC python
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

#+END_SRC

*** Vasp.get_local_potential
args = (self)

Returns x, y, z, and local potential arrays

    We multiply the data by the volume because we are reusing the
    charge density code which divides by volume.


Monkey-patch defined in vasp/getters.py at line 349

[[./vasp/getters.py::349]]
#+BEGIN_SRC python
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

#+END_SRC

*** Vasp.get_neb
args = (self, npi=1)

Returns images, energies if available or runs the job.

    npi = cores per image for running the calculations. Default=1

    show: if True show an NEB plot


Monkey-patch defined in vasp/neb.py at line 46

[[./vasp/neb.py::46]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
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
        # this creates the directories and files if needed.  write out
        # all the images, including initial and final
        if not os.path.isdir(self.directory):
            os.makedirs(self.directory)

        self.set(images=len(self.neb) - 2)
        self.write_incar()
        self.write_kpoints()
        self.write_potcar()
        self.write_db()

        for i, atoms in enumerate(self.neb):
            # zero-padded directory name
            image_dir = os.path.join(self.directory, str(i).zfill(2))
            if not os.path.isdir(image_dir):
                # create if needed.
                os.makedirs(image_dir)
                write_vasp('{0}/POSCAR'.format(image_dir),
                           atoms[self.resort],
                           symbol_count=self.symbol_count)

        # The first and last images need to have real calculators on
        # them so we can write out a DB entry. We need this so we can
        # get the energies on the end-points. Otherwise, there doesn't
        # seem to be a way to do that short of cloning the whole
        # calculation into the end-point directories.

        self.write_db(self.neb[0],
                      os.path.join(self.directory,
                                   '00/DB.db'))

        self.write_db(self.neb[-1],
                      os.path.join(self.directory,
                                   '0{}/DB.db'.format(len(self.neb) - 1)))

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

#+END_SRC

*** Vasp.get_number_of_spins
args = (self)

Returns number of spins.
    1 if not spin-polarized
    2 if spin-polarized



Monkey-patch defined in vasp/getters.py at line 131

[[./vasp/getters.py::131]]
#+BEGIN_SRC python
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

#+END_SRC

*** Vasp.get_occupation_numbers
args = (self, kpt=0, spin=0)

Return the occupation of each k-point.

Monkey-patch defined in vasp/getters.py at line 87

[[./vasp/getters.py::87]]
#+BEGIN_SRC python
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

#+END_SRC

*** Vasp.get_pseudopotentials
args = (self)

Return list of (symbol, path, git-hash) for each POTCAR.

Monkey-patch defined in vasp/getters.py at line 483

[[./vasp/getters.py::483]]
#+BEGIN_SRC python
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

#+END_SRC

*** Vasp.get_state
args = (self)

Determine calculation state based on directory contents.

        Returns an integer for the state.



[[./vasp/vasp_core.py::593]]
#+BEGIN_SRC python
    def get_state(self):
        """Determine calculation state based on directory contents.

        Returns an integer for the state.

        """

        base_input = [os.path.exists(os.path.join(self.directory, f))
                      for f in ['INCAR', 'POSCAR', 'POTCAR', 'KPOINTS']]

        # Check for NEB first.
        if (np.array([os.path.exists(os.path.join(self.directory, f))
                      for f in ['INCAR', 'POTCAR', 'KPOINTS']]).all()
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

#+END_SRC

*** Vasp.get_valence_electrons
args = (self)

Return the number of valence electrons for the atoms.
    Calculated from the POTCAR file.


Monkey-patch defined in vasp/getters.py at line 263

[[./vasp/getters.py::263]]
#+BEGIN_SRC python
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

#+END_SRC

*** Vasp.get_vibrational_frequencies
args = (self)

Returns an array of frequencies in wavenumbers.

    You should have run the calculation already. This function does not
    run a calculation.


Monkey-patch defined in vasp/vib.py at line 155

[[./vasp/vib.py::155]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def get_vibrational_frequencies(self):
    """Returns an array of frequencies in wavenumbers.

    You should have run the calculation already. This function does not
    run a calculation.
    """
    self.update()
    atoms = self.get_atoms()
    N = len(atoms)

    frequencies = []

    f = open(os.path.join(self.directory, 'OUTCAR'), 'r')
    while True:
        line = f.readline()
        if line.startswith(' Eigenvectors and eigenvalues'
                           ' of the dynamical matrix'):
            break
    f.readline()  # skip ------
    f.readline()  # skip two blank lines
    f.readline()
    for i in range(3 * N):
        # the next line contains the frequencies
        line = f.readline()
        fields = line.split()

        if 'f/i=' in line:  # imaginary frequency
            # frequency in wave-numbers
            frequencies.append(complex(float(fields[6]), 0j))
        else:
            frequencies.append(float(fields[7]))
        # now skip 1 one line, a line for each atom, and a blank line
        for j in range(1 + N + 1):
            f.readline()  # skip the next few lines
    f.close()
    return frequencies

#+END_SRC

*** Vasp.get_vibrational_modes
args = (self, mode=None, massweighted=False, show=False, npoints=30, amplitude=0.5)

Read the OUTCAR and get the eigenvectors. Return value depends
    on the arguments.

    mode= None returns all modes
    mode= 2 returns mode 2
    mode=[1, 2] returns modes 1 and 2

    massweighted = True returns sqrt(mass) weighted
    eigenvectors. E.g. M * evectors * M

    show=True makes a trajectory that can be visualized
    npoints = number of points in the trajectory
    amplitude = magnitude of the vibrations

    some special cases to handle:
    ibrion=5 + selective dynamics
       may lead to unexpected number of modes

    if nwrite=3, there will be a sqrt(mass) weighted vectors
    and two sets of vectors.

    I am not sure if these eigenvectors are mass-weighted. And I am
    not sure if the order of the eigenvectors in OUTCAR is the same as
    the atoms.

    Note: it seems like it might be much easier to get this out of
    vasprun.xml


Monkey-patch defined in vasp/vib.py at line 13

[[./vasp/vib.py::13]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def get_vibrational_modes(self,
                          mode=None,
                          massweighted=False,
                          show=False,
                          npoints=30,
                          amplitude=0.5):

    """Read the OUTCAR and get the eigenvectors. Return value depends
    on the arguments.

    mode= None returns all modes
    mode= 2 returns mode 2
    mode=[1, 2] returns modes 1 and 2

    massweighted = True returns sqrt(mass) weighted
    eigenvectors. E.g. M * evectors * M

    show=True makes a trajectory that can be visualized
    npoints = number of points in the trajectory
    amplitude = magnitude of the vibrations

    some special cases to handle:
    ibrion=5 + selective dynamics
       may lead to unexpected number of modes

    if nwrite=3, there will be a sqrt(mass) weighted vectors
    and two sets of vectors.

    I am not sure if these eigenvectors are mass-weighted. And I am
    not sure if the order of the eigenvectors in OUTCAR is the same as
    the atoms.

    Note: it seems like it might be much easier to get this out of
    vasprun.xml
    """
    self.update()

    atoms = self.get_atoms()

    if hasattr(atoms, 'constraints') and self.parameters['ibrion'] == 5:
        # count how many modes to get.
        NMODES = 0
        f = open(os.path.join(self.directory, 'OUTCAR'))
        for line in f:
            if ('f' in line and 'THz' in line and 'cm-1' in line):
                NMODES += 1
        f.close()
    else:
        NMODES = 3 * len(atoms)

    frequencies, eigenvectors = [], []

    # now we find where the data starts. I think the unweighted
    # vectors always come first. if nwrite=3, then there are
    # sqrt(mass) weighted vectors that follow this section

    f = open(os.path.join(self.directory, 'OUTCAR'), 'r')
    while True:
        line = f.readline()
        if line.startswith(' Eigenvectors and eigenvalues'
                           ' of the dynamical matrix'):
            break
    f.readline()   # skip ------
    f.readline()   # skip two blank lines
    f.readline()

    for i in range(NMODES):
        freqline = f.readline()
        fields = freqline.split()

        if 'f/i=' in freqline:  # imaginary frequency
            frequencies.append(complex(float(fields[-2]) * 0.001, 0j))
        else:
            frequencies.append(float(fields[-2]) * 0.001)
        #        X         Y         Z           dx          dy          dz
        f.readline()
        thismode = []
        for i in range(len(atoms)):
            line = f.readline().strip()
            X, Y, Z, dx, dy, dz = [float(x) for x in line.split()]
            thismode.append(np.array([dx, dy, dz]))
        f.readline()  # blank line

        thismode = np.array(thismode)
        # now we need to resort the vectors in this mode so they match
        # the atoms order
        thismode = thismode[self.resort]

        if massweighted:
            # construct M
            numbers = [a.get('number') for a in atoms]
            M = []
            for i in range(len(atoms)):
                for j in range(3):
                    an = numbers[i]
                    M.append(1. / np.sqrt(atomic_masses[an]))
            M = np.array(M)
            M = np.diag(M)  # diagonal array

            thismode = np.dot(M, thismode.flat)

            thismode = thismode.reshape((len(atoms), 3))
        # renormalize the mode
        mag = np.linalg.norm(thismode)
        thismode /= mag

        eigenvectors.append(thismode)
    f.close()

    eigenvectors = np.array(eigenvectors)

    if mode is None:
        retval = (frequencies, eigenvectors)
    else:
        retval = (frequencies[mode], eigenvectors[mode])

    if show:
        from ase.visualize import view
        if mode is None:
            mode = [0]
        elif not isinstance(mode, list):
            mode = [mode]  # make a list for next code

        # symmetric path from -1 to 1 to -1
        X = np.append(np.linspace(0, 1, npoints / 3),
                      np.linspace(1, -1, npoints / 3))
        X = np.append(X,
                      np.linspace(-1, 0, npoints / 3))
        X *= amplitude

        for m in mode:
            traj = []
            for i, x in enumerate(X):
                a = atoms.copy()
                a.positions += x * eigenvectors[m]
                traj += [a]

            view(traj)
    return retval

#+END_SRC

*** Vasp.get_volumetric_data
args = (self, filename=None, **kwargs)

Read filename to read the volumetric data in it.
    Supported filenames are CHG, CHGCAR, and LOCPOT.


Monkey-patch defined in vasp/getters.py at line 282

[[./vasp/getters.py::282]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def get_volumetric_data(self, filename=None, **kwargs):
    """Read filename to read the volumetric data in it.
    Supported filenames are CHG, CHGCAR, and LOCPOT.
    """
    self.update()
    if filename is None:
        filename = os.path.join(self.directory, 'CHG')

    from VaspChargeDensity import VaspChargeDensity

    atoms = self.get_atoms()
    vd = VaspChargeDensity(filename)

    data = np.array(vd.chg)
    n0, n1, n2 = data[0].shape

    # This is the old code, but it doesn't seem to work anymore.
    # s0 = np.linspace(0, 1, num=n0, endpoint=False)
    # s1 = np.linspace(0, 1, num=n1, endpoint=False)
    # s2 = np.linspace(0, 1, num=n2, endpoint=False)

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

#+END_SRC

*** Vasp.in_queue
args = (self)

Return True or False if the directory has a job in the queue.

Monkey-patch defined in vasp/runner.py at line 30

[[./vasp/runner.py::30]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def in_queue(self):
    """Return True or False if the directory has a job in the queue."""
    if self.get_db('jobid') is None:
        log.debug('jobid not found for calculation.')
        return False
    else:
        # get the jobid
        jobid = self.get_db('jobid')
        # see if jobid is in queue
        _, jobids_in_queue, _ = getstatusoutput('qselect',
                                                stdout=subprocess.PIPE,
                                                stderr=subprocess.PIPE)

        if str(jobid) in jobids_in_queue.split('\n'):
            # get details on specific jobid in case it is complete
            status, output, err = getstatusoutput(['qstat', jobid],
                                                  stdout=subprocess.PIPE,
                                                  stderr=subprocess.PIPE)
            if status == 0:
                lines = output.split('\n')
                fields = lines[2].split()
                job_status = fields[4]
                if job_status == 'C':
                    return False
                else:
                    return True
        else:
            return False

#+END_SRC

*** Vasp.plot_neb
args = (self, show=True)

Return a list of the energies and atoms objects for each image in

    the band.

    by default shows the plot figure


Monkey-patch defined in vasp/neb.py at line 159

[[./vasp/neb.py::159]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
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
    f = interp1d(range(len(energies)),
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
    s = ['$\Delta E$ = {0:1.3f} eV'.format(float(energies[-1]
                                                 - energies[0])),
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

#+END_SRC

*** Vasp.read
args = (self, restart=None)

Read the files in a calculation if they exist.

    restart is ignored, but part of the signature for ase. I am not
    sure what we could use it for.

    sets self.parameters and atoms.



Monkey-patch defined in vasp/readers.py at line 220

[[./vasp/readers.py::220]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def read(self, restart=None):
    """Read the files in a calculation if they exist.

    restart is ignored, but part of the signature for ase. I am not
    sure what we could use it for.

    sets self.parameters and atoms.

    """

    self.neb = None
    # NEB is special and handled separately
    if self.get_state() == vasp.Vasp.NEB:
        self.read_neb()
        return

    # Else read a regular calculation. we start with reading stuff
    # that is independent of the calculation state.
    self.parameters = Parameters()

    if os.path.exists(self.incar):
        self.parameters.update(self.read_incar())
    if os.path.exists(self.potcar):
        self.parameters.update(self.read_potcar())
    if os.path.exists(self.kpoints):
        self.parameters.update(self.read_kpoints())

    # We have to figure out the xc that was used based on the
    # Parameter keys.  We sort the possible xc dictionaries so the
    # ones with the largest number of keys are compared first. This is
    # to avoid false matches of xc's with smaller number of equal
    # keys.
    xc_keys = sorted(vasp.Vasp.xc_defaults,
                     key=lambda k: len(vasp.Vasp.xc_defaults[k]),
                     reverse=True)

    for ex in xc_keys:
        pd = {k: self.parameters.get(k, None)
              for k in vasp.Vasp.xc_defaults[ex]}
        if pd == vasp.Vasp.xc_defaults[ex]:
            self.parameters['xc'] = ex
            break

    # reconstruct ldau_luj. special setups might break this.
    if 'ldauu' in self.parameters:
        ldaul = self.parameters['ldaul']
        ldauj = self.parameters['ldauj']
        ldauu = self.parameters['ldauu']

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

        self.parameters['ldau_luj'] = ldau_luj

    # Now for the atoms. This does depend on the state. self.resort
    # needs to be a list for shuffling constraints if they exist.
    self.resort = self.get_db('resort')
    if self.resort is not None:
        self.resort = list(self.resort)

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
        atoms = ase.io.read(contcar)
    elif os.path.exists(poscar):
        atoms = ase.io.read(poscar)
    else:
        atoms = None

    if atoms is not None:
        atoms = atoms[self.resort]
        self.sort_atoms(atoms)

    self.read_results()

#+END_SRC

*** Vasp.read_incar
args = (self, fname=None)

Read fname (defaults to INCAR).

    Returns a Parameters dictionary from the INCAR.

    This only reads simple INCAR files, e.g. one tag per line, and
    with no comments in the line. There is no knowledge of any Vasp
    keywords in this, and the values are converted to Python types by
    some simple rules.



Monkey-patch defined in vasp/readers.py at line 28

[[./vasp/readers.py::28]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def read_incar(self, fname=None):
    """Read fname (defaults to INCAR).

    Returns a Parameters dictionary from the INCAR.

    This only reads simple INCAR files, e.g. one tag per line, and
    with no comments in the line. There is no knowledge of any Vasp
    keywords in this, and the values are converted to Python types by
    some simple rules.

    """

    if fname is None:
        fname = self.incar

    params = Parameters()

    with open(fname) as f:
        lines = f.readlines()

    # The first line is a comment
    for line in lines[1:]:
        line = line.strip()
        if ";" in line:
            raise Exception('; found. that is not supported.')
        if '#' in line:
            raise Exception('# found. that is not supported.')
        if line == '':
            continue

        key, val = line.split('=')
        key = key.strip().lower()
        val = val.strip()
        # now we need some logic
        if val == '.TRUE.':
            val = True
        elif val == '.FALSE.':
            val = False
        # Match integers by a regexp that includes signs
        # val.isdigit() does not get negative integers right.
        elif re.match('^[-+]?\d+$', val):
            val = int(val)
        elif isfloat(val):
            val = float(val)
        elif len(val.split(' ')) > 1:
            # this is some kind of list separated by spaces
            val = val.split(' ')
            val = [int(x) if re.match('^[-+]?\d+$', x)
                   else float(x) for x in val]
        else:
            # I guess we have a string here.
            pass

        # make sure magmom is returned as a list. This is an issue
        # when there is only one atom. Then it looks like a float.
        if key == 'magmom':
            if not isinstance(val, list):
                val = [val]

        params[key] = val

    return params

#+END_SRC

*** Vasp.read_kpoints
args = (self, fname=None)

Read KPOINTS file.

    Returns a Parameters object of kpoint tags.



Monkey-patch defined in vasp/readers.py at line 93

[[./vasp/readers.py::93]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def read_kpoints(self, fname=None):
    """Read KPOINTS file.

    Returns a Parameters object of kpoint tags.

    """

    if fname is None:
        fname = self.kpoints

    with open(fname) as f:
        lines = f.readlines()

    params = Parameters()

    # first line is a comment
    # second line is the number of kpoints or 0 for automatic kpoints
    nkpts = int(lines[1].strip())

    # third line you have to specify whether the coordinates are given
    # in cartesian or reciprocal coordinates if nkpts is greater than
    # zero. Only the first character of the third line is
    # significant. The only key characters recognized by VASP are 'C',
    # 'c', 'K' or 'k' for switching to cartesian coordinates, any
    # other character will switch to reciprocal coordinates.
    #
    # if nkpts = 0 then the third line will start with m or g for
    # Monkhorst-Pack and Gamma. if it does not start with m or g, an
    # alternative mode is entered that we do not support yet.

    ktype = lines[2].split()[0].lower()[0]
    if nkpts <= 0:
        # automatic mode
        if ktype not in ['g', 'm']:
            raise NotImplementedError('Only Monkhorst-Pack and '
                                      'gamma centered grid supported '
                                      'for restart.')
        if ktype == 'g':
            line5 = np.array([float(lines[4].split()[i]) for i in range(3)])
            if (line5 == np.array([0.0, 0.0, 0.0])).all():
                params['gamma'] = True
            else:
                params['gamma'] = line5

        kpts = [int(lines[3].split()[i]) for i in range(3)]
        params['kpts'] = kpts
    elif nkpts > 0:
        # list of kpts provided. Technically c,k are supported and
        # anything else means reciprocal coordinates.
        if ktype in ['c', 'k', 'r']:
            kpts = []
            for i in range(3, 3 + nkpts):
                # the kpts also have a weight attached to them
                kpts.append([float(lines[i].split()[j])
                             for j in range(4)])
            params['kpts'] = kpts
        # you may also be in line-mode
        elif ktype in ['l']:
            if lines[3][0].lower() == 'r':
                params['reciprocal'] = True

            params['kpts_nintersections'] = nkpts

            kpts = []
            for i in range(4, len(lines)):
                if lines[i] == '':
                    continue
                else:
                    kpts.append([float(lines[i].split()[j])
                                 for j in range(3)])
        else:
            raise NotImplementedError('ktype = %s' % lines[2])

    if ktype == 'r':
        params['reciprocal'] = True

    params['kpts'] = kpts
    return params

#+END_SRC

*** Vasp.read_neb
args = (self)

Read an NEB calculator.

Monkey-patch defined in vasp/readers.py at line 383

[[./vasp/readers.py::383]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def read_neb(self):
    """Read an NEB calculator."""
    import ase
    import glob
    atoms = []
    atoms += [ase.io.read('{}/00/POSCAR'.format(self.directory))]
    for p in glob.glob('{}/0[0-9]/CONTCAR'.format(self.directory)):
        atoms += [ase.io.read(p)]
    atoms += [ase.io.read('{}/0{}/POSCAR'.format(self.directory,
                                                 len(atoms)))]
    self.neb = atoms
    self.parameters = {}
    self.set(images=(len(atoms) - 2))
    self.atoms = atoms[0].copy()

    if os.path.exists(self.incar):
        self.parameters.update(self.read_incar())
    if os.path.exists(self.potcar):
        self.parameters.update(self.read_potcar())
    if os.path.exists(self.kpoints):
        self.parameters.update(self.read_kpoints())

    # Update the xc functional
    xc_keys = sorted(vasp.Vasp.xc_defaults,
                     key=lambda k: len(vasp.Vasp.xc_defaults[k]),
                     reverse=True)

    for ex in xc_keys:
        pd = {k: self.parameters.get(k, None)
              for k in vasp.Vasp.xc_defaults[ex]}
        if pd == vasp.Vasp.xc_defaults[ex]:
            self.parameters['xc'] = ex
            break

#+END_SRC

*** Vasp.read_potcar
args = (self, fname=None)

Read the POTCAR file to get the pp and setups.

    Returns a Parameters dictionary of pp and setups.



Monkey-patch defined in vasp/readers.py at line 174

[[./vasp/readers.py::174]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def read_potcar(self, fname=None):
    """Read the POTCAR file to get the pp and setups.

    Returns a Parameters dictionary of pp and setups.

    """

    if fname is None:
        fname = self.potcar

    params = Parameters()

    potcars = []
    with open(fname) as f:
        lines = f.readlines()

    # first potcar
    potcars += [lines[0].strip()]

    for i, line in enumerate(lines):
        if 'LEXCH  = PE' in line:
            params['pp'] = 'PBE'
        elif 'LEXCH  = CA' in line:
            params['pp'] = 'LDA'
        elif 'LEXCH  = 91' in line:
            params['pp'] = 'GGA'

        if 'End of Dataset' in line and i != len(lines) - 1:
            potcars += [lines[i + 1].strip()]

    potcars = [(x[0], x[1], x[2]) for x in
               [potcar.split() for potcar in potcars]]

    special_setups = []
    for xc, sym, date in potcars:
        if '_' in sym:  # we have a special setup
            symbol, setup = sym.split('_')
            special_setups += [[symbol, '_' + setup]]

    if special_setups:
        params['setups'] = special_setups

    return params

#+END_SRC

*** Vasp.read_results
args = (self)

Read energy, forces, stress, magmom and magmoms from output file.

    Other quantities will be read by other functions. This depends on
    state.



Monkey-patch defined in vasp/readers.py at line 316

[[./vasp/readers.py::316]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def read_results(self):
    """Read energy, forces, stress, magmom and magmoms from output file.

    Other quantities will be read by other functions. This depends on
    state.

    """
    state = self.get_state()
    if state == vasp.Vasp.NEB:
        # This is handled in self.read()
        return

    if state != vasp.Vasp.FINISHED:
        self.results = {}
    else:
        # regular calculation that is finished
        from ase.io.vasp import read_vasp_xml
        if not os.path.exists(os.path.join(self.directory,
                                           'vasprun.xml')):
            exc = 'No vasprun.xml in {}'.format(self.directory)
            raise exceptions.VaspNotFinished(exc)

        atoms = read_vasp_xml(os.path.join(self.directory,
                                           'vasprun.xml')).next()

        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()  # needs to be resorted
        stress = atoms.get_stress()

        resort = self.get_db('resort')
        if self.atoms is None:
            atoms = atoms[resort]
            self.sort_atoms(atoms)
            self.atoms.set_calculator(self)
        else:
            # update the atoms
            self.atoms.positions = atoms.positions[resort]
            self.atoms.cell = atoms.cell
            imm = self.parameters.get('magmom',
                                      [0 for atom in self.atoms])
            self.atoms.set_initial_magnetic_moments(imm)

        self.results['energy'] = energy
        self.results['forces'] = forces[self.resort]
        self.results['stress'] = stress
        self.results['dipole'] = None
        self.results['charges'] = np.array([None for atom in self.atoms])

        magnetic_moment = 0
        magnetic_moments = np.zeros(len(atoms))
        if self.parameters.get('ispin', 0) == 2:
            lines = open(os.path.join(self.directory, 'OUTCAR'),
                         'r').readlines()
            for n, line in enumerate(lines):
                if line.rfind('number of electron  ') > -1:
                    magnetic_moment = float(line.split()[-1])

                if line.rfind('magnetization (x)') > -1:
                    for m in range(len(atoms)):
                        val = float(lines[n + m + 4].split()[4])
                        magnetic_moments[m] = val

        self.results['magmom'] = magnetic_moment
        self.results['magmoms'] = np.array(magnetic_moments)[self.resort]

#+END_SRC

*** Vasp.reset
args = (self)

overwrite to avoid killing self.atoms.

[[./vasp/vasp_core.py::468]]
#+BEGIN_SRC python
    def reset(self):
        """overwrite to avoid killing self.atoms."""
        self.results = {}

#+END_SRC

*** Vasp.run
args = (self)

Convenience function to run calculation.

[[./vasp/vasp_core.py::563]]
#+BEGIN_SRC python
    def run(self):
        """Convenience function to run calculation."""
        return self.potential_energy

#+END_SRC

*** Vasp.set
args = (self, **kwargs)

Set parameters with keyword=value pairs.

    calc.set(xc='PBE')

    A few special kwargs are handled separately to expand them
    prior to setting the parameters. This is done to enable one
    set to track changes.



Monkey-patch defined in vasp/setters.py at line 17

[[./vasp/setters.py::17]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
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

    if 'nsw' in kwargs:
        kwargs.update(self.set_nsw_dict(kwargs['nsw']))

    changed_parameters = FileIOCalculator.set(self, **kwargs)
    if changed_parameters:
        self.reset()
    return changed_parameters

#+END_SRC

*** Vasp.set_ispin_dict
args = (self, val)

Returns dictionary of changes for ispin change.

Monkey-patch defined in vasp/setters.py at line 47

[[./vasp/setters.py::47]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def set_ispin_dict(self, val):
    """Returns dictionary of changes for ispin change."""
    # there are two ways to get magmom in.
    # 1. if you use magmom as a keyword, they are used.
    # 2. if you set magmom on each atom in an Atoms object and do not use
    # magmom then we use the atoms magmom, if we have ispin=2 set.
    # we set lorbit to 11 if ispin=2 so we can get the individual moments.
    if val is None:
        d = {}
        for key in ['ispin', 'magmom', 'lorbit']:
            if key in self.parameters:
                d[key] = None
        return d
    elif val == 1:
        d = {'ispin': 1}
        if 'magmom' in self.parameters:
            d['magmom'] = None
        if 'lorbit' in self.parameters:
            d['lorbit'] = None
        return d
    elif val == 2:
        d = {'ispin': 2}
        if 'magmom' not in self.parameters:
            d['magmom'] = [atom.magmom for atom
                            in self.atoms[self.resort]]
        # set individual magnetic moments.
        if 'lorbit' not in self.parameters:
            d['lorbit'] = 11

        return d

#+END_SRC

*** Vasp.set_label
args = (self, label)

Set working directory.

        In VASP there is no prefix, only the working directory.



[[./vasp/vasp_core.py::383]]
#+BEGIN_SRC python
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

#+END_SRC

*** Vasp.set_ldau_luj_dict
args = (self, val)

Set the ldau_luj parameters.

Monkey-patch defined in vasp/setters.py at line 95

[[./vasp/setters.py::95]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def set_ldau_luj_dict(self, val):
    """Set the ldau_luj parameters."""
    if 'setups' in self.parameters:
        raise Exception('setups and ldau_luj is not supported.')

    if not hasattr(self, 'ppp_list'):
        atoms = self.get_atoms()
        self.sort_atoms(atoms)

    if val is not None:
        atom_types = [x[0] if isinstance(x[0], str)
                      else self.atoms[x[0]].symbol
                      for x in self.ppp_list]

        d = {}

        d['ldaul'] = [val[sym]['L'] for sym in atom_types]
        d['ldauu'] = [val[sym]['U'] for sym in atom_types]
        d['ldauj'] = [val[sym]['J'] for sym in atom_types]
        return d
    else:
        d = {}
        d['ldaul'] = None
        d['ldauu'] = None
        d['ldauj'] = None
        return d

#+END_SRC

*** Vasp.set_nbands
args = (self, N=None, f=1.5)

Convenience function to set NBANDS to N or automatically compute
    nbands for non-spin-polarized calculations.

    nbands = int(nelectrons/2 + nions*f)

    this formula is suggested at
    http://cms.mpi.univie.ac.at/vasp/vasp/NBANDS_tag.html for
    transition metals f may be as high as 2.



Monkey-patch defined in vasp/setters.py at line 124

[[./vasp/setters.py::124]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def set_nbands(self, N=None, f=1.5):
    """Convenience function to set NBANDS to N or automatically compute
    nbands for non-spin-polarized calculations.

    nbands = int(nelectrons/2 + nions*f)

    this formula is suggested at
    http://cms.mpi.univie.ac.at/vasp/vasp/NBANDS_tag.html for
    transition metals f may be as high as 2.

    """
    if N is not None:
        self.set(nbands=int(N))
        return
    atoms = self.get_atoms()
    nelectrons = self.get_valence_electrons()
    nbands = int(np.ceil(nelectrons / 2.) + len(atoms) * f)
    self.set(nbands=nbands)

#+END_SRC

*** Vasp.set_nsw_dict
args = (self, val)

Set nsw parameter.

    The default lwave behavior is False, but if nsw > 0 it makes sense
    to turn it on in case of restarts.



Monkey-patch defined in vasp/setters.py at line 145

[[./vasp/setters.py::145]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def set_nsw_dict(self, val):
    """Set nsw parameter.

    The default lwave behavior is False, but if nsw > 0 it makes sense
    to turn it on in case of restarts.

    """

    d = {'nsw': val}

    if val > 0:
        d['lwave'] = True
    elif val == 0:
        d['lwave'] = False
    else:
        d['lwave'] = False
    return d

#+END_SRC

*** Vasp.set_rwigs_dict
args = (self, val)

Return rwigs parameters.

Monkey-patch defined in vasp/setters.py at line 80

[[./vasp/setters.py::80]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def set_rwigs_dict(self, val):
    """Return rwigs parameters."""
    d = {}
    if val is None:
        d['rwigs'] = None
        d['lorbit'] = None
    else:
        # val is a dictionary {sym: rwigs}
        # rwigs needs to be in the order of the potcars
        d['rwigs'] = [val[x[0]] for x in self.ppp_list]

    return d

#+END_SRC

*** Vasp.set_xc_dict
args = (self, val)

Set xc parameter.

    Adds all the xc_defaults flags for the chosen xc.



Monkey-patch defined in vasp/setters.py at line 165

[[./vasp/setters.py::165]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def set_xc_dict(self, val):
    """Set xc parameter.

    Adds all the xc_defaults flags for the chosen xc.

    """
    d = {'xc': val.lower()}
    oxc = self.parameters.get('xc', None)
    if oxc:
        for key in vasp.Vasp.xc_defaults[oxc.lower()]:
            if key in self.parameters:
                d[key] = None
    d.update(vasp.Vasp.xc_defaults[val.lower()])
    return d

#+END_SRC

*** Vasp.sort_atoms
args = (self, atoms=None)

Generate resort list, and make list of POTCARs to use.

        Returns None.



[[./vasp/vasp_core.py::282]]
#+BEGIN_SRC python
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

        return atoms[self.resort]

#+END_SRC

*** Vasp.stop_if
args = (self, condition=None)

Stop program if condition is truthy.

[[./vasp/vasp_core.py::553]]
#+BEGIN_SRC python
    def stop_if(self, condition=None):
        """Stop program if condition is truthy."""
        if condition:
            import sys
            sys.exit()

#+END_SRC

*** Vasp.update
args = (self, atoms=None)

Updates calculator.

        If a calculation is required,  run it, otherwise updates results.



[[./vasp/vasp_core.py::472]]
#+BEGIN_SRC python
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

#+END_SRC

*** Vasp.view
args = (self, index=None)

Visualize the calculation.



[[./vasp/vasp_core.py::711]]
#+BEGIN_SRC python
    def view(self, index=None):
        """Visualize the calculation.

        """
        from ase.visualize import view
        if index is not None:
            return view(self.traj[index])
        else:
            return view(self.traj)

#+END_SRC

*** Vasp.wait
args = (self)

Stop program if not ready.

[[./vasp/vasp_core.py::559]]
#+BEGIN_SRC python
    def wait(self):
        """Stop program if not ready."""
        self.stop_if(self.potential_energy is None)

#+END_SRC

*** Vasp.write_db
args = (self, atoms=None, fname=None, data=None, **kwargs)

Write the DB file.

    atoms can be any atoms object, defaults to self.get_atoms().
    fname can be anything, defaults to self.directory/DB.db

    data is a dictionary of data to store.

    kwargs is key=value pairs to store with the atoms.

    Existing data and kwargs are preserved. You can delete kwargs by
    setting them to None. You can delete data by setting the key to
    None in the data dictionary.

    Only row 1 should be in this database.



Monkey-patch defined in vasp/writers.py at line 31

[[./vasp/writers.py::31]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def write_db(self, atoms=None, fname=None, data=None, **kwargs):
    """Write the DB file.

    atoms can be any atoms object, defaults to self.get_atoms().
    fname can be anything, defaults to self.directory/DB.db

    data is a dictionary of data to store.

    kwargs is key=value pairs to store with the atoms.

    Existing data and kwargs are preserved. You can delete kwargs by
    setting them to None. You can delete data by setting the key to
    None in the data dictionary.

    Only row 1 should be in this database.

    """
    from ase.db import connect

    if fname is None:
        fname = os.path.join(self.directory, 'DB.db')

    fdata = {'resort': self.resort,
             'parameters': self.parameters,
             'ppp_list': self.ppp_list}
    fkv = {'path': self.directory}

    # get current data and keywords
    if os.path.exists(fname):

        with connect(fname) as con:
            try:
                temp_atoms = con.get_atoms(id=1,
                                           add_additional_information=True)
                fdata.update(temp_atoms.info['data'])
                fkv.update(temp_atoms.info['key_value_pairs'])
            except KeyError:
                pass

    # Update fdata from input args. None removes keywords and data
    # elements
    if data is not None:
        for key, val in data.iteritems():
            if val is None and key in fdata:
                del fdata[key]
            else:
                fdata.update({key: val})

    # update key-value pairs from input args
    for key, val in kwargs.iteritems():
        # we use None to delete keys
        if val is None and key in fkv:
            del fkv[key]
        elif val is not None:
            fkv.update({key: val})

    if atoms is None:
        atoms = self.get_atoms()

    # TODO: NEB? should contain images?
    # write out in non-append mode.
    with connect(fname, append=False) as con:
        con.write(atoms,
                  data=fdata,
                  **fkv)

#+END_SRC

*** Vasp.write_incar
args = (self, incar=None)

Writes out the INCAR file.

    Boolean values are written as .TRUE./.FALSE.
    integers/floats and strings are written out as is
    lists/tuples are written out as space separated values/



Monkey-patch defined in vasp/writers.py at line 111

[[./vasp/writers.py::111]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def write_incar(self, incar=None):
    """Writes out the INCAR file.

    Boolean values are written as .TRUE./.FALSE.
    integers/floats and strings are written out as is
    lists/tuples are written out as space separated values/

    """

    if incar is None:
        incar = os.path.join(self.directory, 'INCAR')

    incar_keys = list(set(self.parameters) - set(self.special_kwargs))
    d = {key: self.parameters[key] for key in incar_keys}

    with open(incar, 'w') as f:
        f.write('INCAR created by Atomic Simulation Environment\n')
        for key, val in d.iteritems():
            key = ' ' + key.upper()
            if val is None:
                # Do not write out None values
                # It is how we delete tags
                pass
            elif isinstance(val, bool):
                s = '.TRUE.' if val else '.FALSE.'
                f.write('{} = {}\n'.format(key, s))
            elif isinstance(val, list) or isinstance(val, tuple):
                s = ' '.join([str(x) for x in val])
                f.write('{} = {}\n'.format(key, s))
            else:
                f.write('{} = {}\n'.format(key, val))

#+END_SRC

*** Vasp.write_input
args = (self, atoms=None, properties=None, system_changes=None)

Writes all input files required for a calculation.

Monkey-patch defined in vasp/writers.py at line 17

[[./vasp/writers.py::17]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def write_input(self, atoms=None, properties=None, system_changes=None):
    """Writes all input files required for a calculation."""
    # this creates the directory if needed
    FileIOCalculator.write_input(self, atoms, properties, system_changes)

    if 'spring' not in self.parameters:  # do not write if NEB
        self.write_poscar()
    self.write_incar()
    self.write_kpoints()
    self.write_potcar()
    self.write_db()

#+END_SRC

*** Vasp.write_kpoints
args = (self, fname=None)

Write out the KPOINTS file.

    The KPOINTS file format is as follows:

    line 1: a comment
    line 2: number of kpoints
        n <= 0   Automatic kpoint generation
        n > 0    explicit number of kpoints
    line 3: kpt format
        if n > 0:
            C,c,K,k = cartesian coordinates
            anything else = reciprocal coordinates
        if n <= 0
            M,m,G,g for Monkhorst-Pack or Gamma grid
            anything else is a special case
    line 4: if n <= 0, the Monkhorst-Pack grid
        if n > 0, then a line per kpoint
    line 5: if n <=0 it is the gamma shift

    After the kpts may be tetrahedra, but we do now support that for
    now.



Monkey-patch defined in vasp/writers.py at line 145

[[./vasp/writers.py::145]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def write_kpoints(self, fname=None):
    """Write out the KPOINTS file.

    The KPOINTS file format is as follows:

    line 1: a comment
    line 2: number of kpoints
        n <= 0   Automatic kpoint generation
        n > 0    explicit number of kpoints
    line 3: kpt format
        if n > 0:
            C,c,K,k = cartesian coordinates
            anything else = reciprocal coordinates
        if n <= 0
            M,m,G,g for Monkhorst-Pack or Gamma grid
            anything else is a special case
    line 4: if n <= 0, the Monkhorst-Pack grid
        if n > 0, then a line per kpoint
    line 5: if n <=0 it is the gamma shift

    After the kpts may be tetrahedra, but we do now support that for
    now.

    """
    if fname is None:
        fname = os.path.join(self.directory, 'KPOINTS')

    p = self.parameters

    kpts = p.get('kpts', None)  # this is a list, or None

    if kpts is None:
        NKPTS = None
    elif len(np.array(kpts).shape) == 1:
        NKPTS = 0  # automatic
    else:
        NKPTS = len(p['kpts'])

    # figure out the mode
    if NKPTS == 0 and not p.get('gamma', None):
        MODE = 'm'  # automatic monkhorst-pack
    elif NKPTS == 0 and p.get('gamma', None):
        MODE = 'g'  # automatic gamma monkhorst pack
    # we did not trigger automatic kpoints
    elif p.get('kpts_nintersections', None) is not None:
        MODE = 'l'
    elif p.get('reciprocal', None) is True:
        MODE = 'r'
    else:
        MODE = 'c'

    with open(fname, 'w') as f:
        # line 1 - comment
        f.write('KPOINTS created by Atomic Simulation Environment\n')
        # line 2 - number of kpts
        if MODE in ['c', 'k', 'm', 'g', 'r']:
            f.write('{}\n'.format(NKPTS))
        elif MODE in ['l']:  # line mode, default intersections is 10
            f.write('{}\n'.format(p.get('kpts_nintersections')))

        # line 3
        if MODE in ['m', 'g']:
            if MODE == 'm':
                f.write('Monkhorst-Pack\n')  # line 3
            elif MODE == 'g':
                f.write('Gamma\n')
        elif MODE in ['c', 'k']:
            f.write('Cartesian\n')
        elif MODE in ['l']:
            f.write('Line-mode\n')
        else:
            f.write('Reciprocal\n')

        # line 4
        if MODE in ['m', 'g']:
            f.write('{0} {1} {2}\n'.format(*p.get('kpts', (1, 1, 1))))
        elif MODE in ['c', 'k', 'r']:
            for n in range(NKPTS):
                # I assume you know to provide the weights
                f.write('{0} {1} {2} {3}\n'.format(*p['kpts'][n]))
        elif MODE in ['l']:
            if p.get('reciprocal', None) is False:
                f.write('Cartesian\n')
            else:
                f.write('Reciprocal\n')
            for n in range(NKPTS):
                f.write('{0} {1} {2}\n'.format(*p['kpts'][n]))

        # line 5 - only if we are in automatic mode
        if MODE in ['m', 'g']:
            if p.get('gamma', None) and len(p.get('gamma')) == 3:
                f.write('{0} {1} {2}\n'.format(*p['gamma']))
            else:
                f.write('0.0 0.0 0.0\n')

#+END_SRC

*** Vasp.write_poscar
args = (self, fname=None)

Write the POSCAR file.

Monkey-patch defined in vasp/writers.py at line 99

[[./vasp/writers.py::99]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def write_poscar(self, fname=None):
    """Write the POSCAR file."""
    if fname is None:
        fname = os.path.join(self.directory, 'POSCAR')

    from ase.io.vasp import write_vasp
    write_vasp(fname,
               self.atoms_sorted,
               symbol_count=self.symbol_count)

#+END_SRC

*** Vasp.write_potcar
args = (self, fname=None)

Writes the POTCAR file.

    POTCARs are expected in $VASP_PP_PATH.



Monkey-patch defined in vasp/writers.py at line 242

[[./vasp/writers.py::242]]
#+BEGIN_SRC python
@monkeypatch_class(vasp.Vasp)
def write_potcar(self, fname=None):
    """Writes the POTCAR file.

    POTCARs are expected in $VASP_PP_PATH.

    """
    if fname is None:
        fname = os.path.join(self.directory, 'POTCAR')

    with open(fname, 'wb') as potfile:
        for _, pfile, _ in self.ppp_list:
            pfile = os.path.join(os.environ['VASP_PP_PATH'], pfile)
            with open(pfile) as f:
                potfile.write(f.read())

#+END_SRC
