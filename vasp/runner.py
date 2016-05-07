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

from ase.calculators.calculator import Calculator,\
    FileIOCalculator


def getstatusoutput(*args, **kwargs):
    """Helper function to replace the old commands.getstatusoutput.

    Returns the returncode, stdout and sterr associated with the command.

    getstatusoutput([command], stdin=subprocess.PIPE)

    """
    p = subprocess.Popen(*args, **kwargs)
    stdout, stderr = p.communicate()
    return (p.returncode, stdout, stderr)


@monkeypatch_class(vasp.Vasp)
def in_queue(self):
    """Return True or False if the directory has a job in the queue."""
    if 'jobid' not in self.metadata:
        log.debug('jobid not found in {}', self.metadata)
        return False
    else:
        # get the jobid
        jobid = self.metadata['jobid']
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
        raise VaspQueued(self.metadata['jobid'])

    # not in queue. Delete the jobid
    if 'jobid' in self.metadata:
        del self.metadata['jobid']
        self.write_metadata()

        # we should check for errors here.
        self.read_results()
        return

    # The subclass implementation should first call this
    # implementation to set the atoms attribute.
    Calculator.calculate(self, atoms, properties, system_changes)

    self.write_input(atoms, properties, system_changes)

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
                    log.debug('queue.nodes = {0}'.format(VASPRC['queue.nodes']))
                    log.debug('queue.ppn = {0}'.format(VASPRC['queue.ppn']))
                    log.debug('multiprocessing.cores_per_process'
                              '= {0}'.format(VASPRC['multiprocessing.cores_per_process']))
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
            vaspcmd = VASPRC['vasp.executable.serial']
            exitcode = os.system(vaspcmd)
            return exitcode
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

    self.metadata['jobid'] = out.strip()
    self.write_metadata()

    raise VaspSubmitted(out)
