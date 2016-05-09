"""Setter functions for the Vasp calculator.

These are special keyword setters that are needed because a special
tag represents some set of keywords, e.g. the xc tag, or a change
implies multiple changes.

The general strategy is to construct a dictionary of the changes that
is passed back to the set function.

"""
import vasp
from monkeypatch import monkeypatch_class


@monkeypatch_class(vasp.Vasp)
def set_ispin_dict(self, val):
    """Returns dictionary of changes for ispin change."""
    # there are two ways to get magmoms in.
    # 1. if you use magmoms as a keyword, they are used.
    # 2. if you set magmom on each atom in an Atoms object and do not use
    # magmoms then we use the atoms magmoms, if we have ispin=2 set.
    if val is None:
        d = {}
        for key in ['ispin', 'magmoms', 'lorbit']:
            if key in self.parameters:
                d[key] = None
            return d
    elif val == 1:
        d = {'ispin': 1}
        if 'magmoms' in self.parameters:
            d['magmoms'] = None
        if 'lorbit' in self.parameters:
            d['lorbit'] = None
        return d
    elif val == 2:
        d = {'ispin': 2}
        if 'magmoms' not in self.parameters:
            d['magmoms'] = [atom.magmom for atom
                            in self.atoms[self.resort]]
        # set individual magnetic moments.
        if 'lorbit' not in self.parameters:
            d['lorbit'] = 11

        return d


@monkeypatch_class(vasp.Vasp)
def set_xc_dict(self, val):
    """Set xc parameter."""
    d = {'xc': val}
    oxc = self.parameters.get('xc', None)
    if oxc:
        for key in vasp.Vasp.xc_defaults[oxc.lower()]:
            if key in self.parameters:
                d[key] = None
    d.update(vasp.Vasp.xc_defaults[val.lower()])
    return d


@monkeypatch_class(vasp.Vasp)
def set_ldau_luj_dict(self, val):
    """Set the ldau_luj parameters."""
    if 'setups' in self.parameters:
        raise Exception('setups and ldau_luj is not supported.')

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
