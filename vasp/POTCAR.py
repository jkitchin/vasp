import re
from subprocess import Popen, PIPE


def get_ZVAL(potcar):
    """Return the ZVAL for a potcar file.

    parse this line:
       POMASS =  106.420; ZVAL   =   10.000    mass and valenz
    """
    # First check if it is a .Z type file
    if potcar.endswith('.Z'):
        cmdlist = ['zcat', potcar]
        p = Popen(cmdlist, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        if out == '' or err != '':
            raise Exception('Cannot read POTCAR.Z:\n\n{0}'.format(err))

        lines = out.split('\n')

    else:
        with open(potcar, 'r') as f:
            lines = f.readlines()

    for line in lines:
        if 'ZVAL' in line:
            m = re.search('ZVAL   =\s*([0-9]*\.?[0-9]*)', line)
            return float(m.group(1))

    return


def get_ENMAX(potcar):
    """ Return ENMAX from the potcar file."""
    with open(potcar) as f:
        for line in f:
            if 'ENMAX' in line:
                m = re.search('ENMAX\s*=\s*(?P<ENMAX>[0-9]+.[0-9]+);', line)
                return float(m.groupdict()['ENMAX'])


def get_ENMIN(potcar):
    """ Return ENMIN from the potcar file."""
    with open(potcar) as f:
        for line in f:
            if 'ENMIN' in line:
                regex = 'ENMIN\s*=\s*(?P<ENMIN>[0-9]+.[0-9]+)\s+eV'
                m = re.search(regex, line)
                return float(m.groupdict()['ENMIN'])
