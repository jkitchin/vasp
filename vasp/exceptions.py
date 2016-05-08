"""Exceptions for Vasp."""


class VaspSubmitted(Exception):
    def __init__(self, jobid):
        self.jobid = jobid

    def __str__(self):
        return str(self.jobid)


class VaspQueued(Exception):
    def __init__(self, msg='Queued', cwd=None):
        self.msg = msg
        self.cwd = cwd

    def __str__(self):
        return '{}'.format(self.msg)


class VaspRunning(Exception):
    pass


class VaspNotFinished(Exception):
    def __init__(self, message=''):
        self.message = message

    def __str__(self):
        return self.message


class VaspEmptyCONTCAR(Exception):
    def __init__(self, message=''):
        self.message = message

    def __str__(self):
        return self.message


class VaspNotConverged(Exception):
    pass


class VaspUnknownState(Exception):
    pass


class VaspWarning(Exception):
    '''Exception class for Vasp warnings that cause problems in Vasp.'''
    def __init__(self, message=''):
        self.message = message

    def __str__(self):
        return self.message
