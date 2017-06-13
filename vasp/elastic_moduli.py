"""Module to get elastic moduli from Vasp calculations."""
import os
from .vasp import Vasp
import numpy as np
from .monkeypatch import monkeypatch_class


@monkeypatch_class(Vasp)
def get_elastic_moduli(self):
    """Returns the total elastic moduli in GPa.

    (i.e. the rigid ion and contributions from relaxation) from the
    OUTCAR file.

    you must run with IBRION=6 and ISIF>= 3 for this output to exist.

    There are also contributions from ionic relaxation
    ELASTIC MODULI CONTR FROM IONIC RELAXATION (kBar)
    and the rigid moduli
    SYMMETRIZED ELASTIC MODULI (kBar)

    For now these are not returned.
    """
    assert self.parameters.get('ibrion', 0) == 6
    assert self.parameters.get('isif', 0) >= 3

    self.update()

    with open(os.path.join(self.directory, 'OUTCAR')) as f:
        lines = f.readlines()

    TEM = []
    for i, line in enumerate(lines):
        if line.startswith(' TOTAL ELASTIC MODULI (kBar)'):
            j = i + 3
            data = lines[j: j + 6]
            break

    for line in data:
        # each line looks like this:
        # XX 2803.5081 1622.6085 1622.6085 0.0000 0.0000 0.0000
        TEM += [[float(x) for x in line.split()[1:]]]

    return np.array(TEM) * 0.1  # (convert kbar to GPa)
