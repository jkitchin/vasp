"""Run Vasp jobs in a queue.

This assumes you use Torque.
"""

import os
import subprocess
import vasp
from vasprc import VASPRC
from vasp import log
from exceptions import VaspSubmitted, VaspQueued
from monkeypatch import monkeypatch_class

from ase.calculators.calculator import Calculator


def getstatusoutput(*args, **kwargs):
    """Helper function to replace the old commands.getstatusoutput.

    Returns the returncode, stdout and sterr associated with the command.

    getstatusoutput([command], stdin=subprocess.PIPE)

    """
    p = subprocess.Popen(*args, **kwargs)
    stdout, stderr = p.communicate()
    return (p.returncode, stdout, stderr)


@monkeypatch_class(vasp.Vasp)
def jobid(self):
    """Return jobid for the calculation."""
    return self.get_db('jobid')


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

    if VASPRC['mode'] is None:
        log.debug('mode is None. not running')
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


@monkeypatch_class(vasp.Vasp)
def set_memory(self,
               eff_loss=0.1):
    """ Sets the recommended memory needed for a VASP calculation

    Code retrieves memory estimate based on the following priority:
    1) DB file
    2) existing OUTCAR
    3) run partial diagnostic calculation

    The final method determines the memory requirements from
    KPOINT calculations run locally before submission to the queue

    returns the memory estimate from the OUTCAR file

    :param eff_loss: Estimated loss in computational efficiency
                     when parallelizing over multiple processors.
                     10% is a safe upper bound from personal experience.
    :type eff_loss: float
    """

    # Attempt to get the recommended memory from DB
    memory = self.get_db('memory')

    if memory is None:
        # Check if an OUTCAR exists from a previous run

        # WARNING: if calculation was run with > 1 ppn, estimate
        # in OUTCAR will reflect that recommended memory for that
        # number of ppn. This can lead to over-estimation of memory
        # requirements on calculations where set_required_memory
        # was not used initially.
        if os.path.exists(os.path.join(self.directory, 'OUTCAR')):
            memory = self.get_memory()

            # Write the recommended memory to the DB file
            self.write_db(data={'memory': memory})

        # If no OUTCAR exists, we run a 'dummy' calculation
        else:
            try:
                original_ialgo = self.parameters['ialgo']
            except(KeyError):
                original_ialgo = None
            self.set(ialgo=-1)

            # Generate the base files needed for VASP calculation
            from ase.calculators.calculator import FileIOCalculator
            FileIOCalculator.write_input(self, None, None, None)
            self.write_poscar()
            self.write_incar()
            if 'kspacing' not in self.parameters:
                self.write_kpoints()
            self.write_potcar()

            # Need to pass a function to Timer for delayed execution
            def kill():
                process.kill()

            # We only need the memory estimate, so we can greatly
            # accelerate the process by terminating after we have it
            cwd = os.getcwd()
            os.chdir(self.directory)
            from subprocess import Popen, PIPE
            process = Popen(VASPRC['vasp.executable.serial'],
                            stdout=PIPE)

            from threading import Timer
            import time
            timer = Timer(20.0, kill)
            timer.start()
            while True:
                if timer.is_alive():
                    memory = self.get_memory()
                    if memory:
                        timer.cancel()
                        process.terminate()
                        break
                    else:
                        time.sleep(0.1)
                else:
                    raise RuntimeError('Memory estimate timed out')

            os.chdir(cwd)
            # return to original settings
            if original_ialgo:
                self.set(ialgo=original_ialgo)
            else:
                del self.parameters['ialgo']
            self.write_incar()

            # Write the recommended memory to the DB file
            self.write_db(data={'memory': memory})

            # Remove all non-initialization files
            files = ['CHG', 'CHGCAR', 'CONTCAR', 'DOSCAR',
                     'EIGENVAL', 'IBZKPT', 'OSZICAR', 'PCDAT',
                     'vasprun.xml', 'OUTCAR', 'WAVECAR', 'XDATCAR']

            for f in files:
                os.unlink(os.path.join(self.directory, f))

    # One node will require the memory read from the OUTCAR
    processors = VASPRC['queue.ppn'] + VASPRC['queue.nodes']

    # Apply eff_loss
    if processors > 1:
        eff_factor = 1 + (eff_loss * float(processors))

    # Rounded up to the nearest GB, and set the memory
    import math
    total_memory = int(math.ceil(eff_factor * float(memory)))
    VASPRC['queue.mem'] = '{0}GB'.format(total_memory)

    # Return the memory as read from the OUTCAR
    return memory


@monkeypatch_class(vasp.Vasp)
def qdel(self, *options):
    """Delete job from the queue.

    options are strings passed to the qdel command.

    This
    >>> calc.qdel('-p')

    is equivalent to the shell-command 'qdel -p jobid'.
    """
    if self.in_queue():
        jobid = self.get_db('jobid')
        cmd = ['qdel'] + list(options) + [jobid]
        status, output, err = getstatusoutput(cmd,
                                              stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE)
        if status != 0:
            print(output + err)
        return status, output
    return '{} not in queue.'.format(self.directory)


@monkeypatch_class(vasp.Vasp)
def qstat(self, *options):
    """Get queue status of the job.

    options are strings, e.g. '-f'

    """
    if self.in_queue():
        jobid = self.get_db('jobid')
        cmd = ['qstat'] + list(options) + [jobid]

        status, output, err = getstatusoutput(cmd,
                                              stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE)
        if status == 0:
            print(output)
        else:
            print(output + err)
    else:
        print('{} not in queue.'.format(self.directory))


@monkeypatch_class(vasp.Vasp)
def qalter(self, *options):
    """Run qalter on the jobid.

    E.g.
    >>> calc.qalter('-l', 'walltime=10:00:00')

    """
    jobid = self.get_db('jobid')
    cmd = ['qalter'] + list(options) + [jobid]

    status, output, err = getstatusoutput(cmd,
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE)
    return status, output


@monkeypatch_class(vasp.Vasp)
def xterm(self):
    """Open an xterm in the calculator directory."""

    cmd = 'xterm -e "cd {}; ls && /bin/bash"'.format(self.directory)
    os.system(cmd)


@monkeypatch_class(vasp.Vasp)
def qoutput(self):
    """Print job output from the queue."""
    jobid = self.jobid()
    ou = os.path.join(self.directory, jobid + '.OU')
    if not self.in_queue() and os.path.exists(ou):
        with open(ou) as f:
            return f.read()
    else:
        return "In queue or no output found."


def torque(cls):
    """Returns an array of jobids and actions suitable for org-mode.

    | directory | jobid | qdel |

    The directory link opens in xterm or dired,
    The jobid runs qstat on the job_status
    qdel will delete the job from the queue.
    """
    jobids = [calc.jobid() for calc in vasp.Vasp.calculators]

    qstat = ['[[shell:qstat {}][{}]]'.format(jobid, jobid)
             for jobid in jobids]
    qdel = ['[[shell:qdel {}][qdel]]'.format(jobid)
            for jobid in jobids]

    dirs = [calc.directory
            for calc in vasp.Vasp.calculators]

    s = '[[shell:xterm -e "cd {}; ls && /bin/bash"][{}]]'
    xterm = [s.format(d, os.path.relpath(d))
             for d in dirs]

    s = '[[elisp:(find-file "{}")][dired]]'
    dired = [s.format(d)
             for d in dirs]

    return '\n'.join(['| {0} {1} | {2} | {3} |'.format(xt, dd, qs, qd)
                      for xt, qs, qd, dd in zip(xterm, qstat, qdel, dired)])

vasp.Vasp.torque = classmethod(torque)
