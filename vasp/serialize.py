"""Properties for serializing vasp calculations."""
import os
from os.path import join
import vasp as VASP


def vasp(self):
    """Return files as strings.

    except the POTCAR.

    """
    s = 'INCAR'
    s += '\n' + len('INCAR') * "-" + '\n'
    with open(join(self.directory, 'INCAR')) as f:
        s += f.read() + '\n\n'

    s += 'POSCAR\n' + '-' * len('POSCAR') + '\n'
    with open(join(self.directory, 'POSCAR')) as f:
        s += f.read() + '\n\n'

    s += 'KPOINTS\n' + '-' * len('KPOINTS') + '\n'
    with open(join(self.directory, 'KPOINTS')) as f:
        s += f.read() + '\n\n'

    s += 'POTCAR\n' + '-' * len('POTCAR') + '\n'
    s += 'cat ' + ' '.join(['$VASP_PP_PATH/' + x[1] for x in self.ppp_list])
    s += ' > POTCAR'

    return s

setattr(VASP.Vasp, 'vasp', property(vasp))


def vasp_json(self):
    """Return a json representation."""
    json = join(self.directory, 'DB.json')

    self.write_db(fname=json)
    with open(json) as f:
        s = f.read()

    os.unlink(json)
    return s

setattr(VASP.Vasp, 'json', property(vasp_json))


def vasp_jsonpp(self):
    """Return a pretty-printed json representation."""
    import json

    jsonf = join(self.directory, 'DB.json')
    self.write_db(fname=jsonf)
    with open(jsonf) as f:
        d = json.loads(f.read())
    os.unlink(jsonf)
    return json.dumps(d, sort_keys=True, indent=4)

setattr(VASP.Vasp, 'jsonpp', property(vasp_jsonpp))
