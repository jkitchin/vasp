# A Mongo database for ASE calculations

"""This module will be like the ase-db but different in the following ways:

1. Booleans are stored as booleans.
2. There is no numeric id.
3. Tags are stored in an array.
"""

import os
import numpy as np
from collections import OrderedDict
import datetime
import json
import hashlib

from pymongo import MongoClient
from ase import Atoms, Atom
from ase.io.jsonio import encode
from vasp import Vasp
import spglib

def mongo_atoms_doc(atoms):
    """Return a dictionary of an Atoms object."""
    d = OrderedDict(atoms=[{'symbol': atom.symbol,
                            'position': json.loads(encode(atom.position)),
                            'tag': atom.tag,
                            'index': atom.index,
                            'charge': atom.charge,
                            'momentum': json.loads(encode(atom.momentum)),
                            'magmom': atom.magmom}
                           for atom in atoms],
                    cell=atoms.cell,
                    pbc=atoms.pbc,
                    info=atoms.info,
                    constraints=[c.todict() for c in atoms.constraints])

    d['natoms'] = len(atoms)
    cell = atoms.get_cell()
    if cell is not None and np.linalg.det(cell) > 0:
        d['volume'] = atoms.get_volume()

    d['mass'] = sum(atoms.get_masses())

    syms = atoms.get_chemical_symbols()
    d['chemical_symbols'] = list(set(syms))
    d['symbol_counts'] = {sym: syms.count(sym) for sym in syms}
    d['spacegroup'] = spglib.get_spacegroup(atoms)

    return json.loads(encode(d))


class MongoDatabase(MongoClient):

    def __init__(self,
                 host='localhost',
                 port=27017,
                 database='ase',
                 collection='atoms',
                 user=None,
                 password=None):
        """
        user and password are currently unused.
        """
        MongoClient.__init__(self, host, port)

        self.db = self[database]
        self.collection = getattr(self.db, collection)

    def write(self, atoms, **kwargs):
        """
        atoms is an ase.atoms.Atoms object.
        kwargs are key-value pairs that will be written to the database.

        Returns the inserted id.
        """

        d = OrderedDict(atoms=mongo_atoms_doc(atoms))

        # Calculated values
        if atoms.get_calculator() is not None:
            # Need some calculator data
            calc = atoms.get_calculator()
            d['calculator'] = calc.todict()

        d['user'] = os.getenv('USER'),
        # This is a has of what is inserted. You might use it to check
        # for uniqueness of the insert.
        d['inserted-hash'] = hashlib.sha1(json.dumps(d)).hexdigest()

        # Created time.
        d['ctime'] = datetime.datetime.utcnow()
        # Modified time - depends on user to update
        d['mtime'] = datetime.datetime.utcnow()

        d.update(kwargs)

        return self.collection.insert_one(d).inserted_id

    def find(self, *args, **kwargs):
        """Thin wrapper for collection.find().

        """
        return self.collection.find(*args, **kwargs)

    def get_atoms(self, *args, **kwargs):
        """Return an atoms object for each match in filter.

        args and kwargs are passed to the collection.find function
        """

        cursor = self.collection.find(*args, **kwargs)
        for doc in cursor:
            atoms = Atoms([Atom(atom['symbol'],
                                atom['position'],
                                tag=atom['tag'],
                                momentum=atom['momentum'],
                                magmom=atom['magmom'],
                                charge=atom['charge'])
                           for atom in doc['atoms']['atoms']],
                          cell=doc['atoms']['cell'])

            calc_data = doc['calculator']
            pars = calc_data['parameters']
            calc = Vasp(calc_data['path'], **pars)
            atoms.set_calculator(calc)

            # TODO the calculator
            yield atoms
