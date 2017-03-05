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
import pickle
from pymongo import MongoClient
from ase import Atoms, Atom


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

        d = OrderedDict(user=os.getenv('USER'),
                        ctime=datetime.datetime.utcnow(),
                        mtime=datetime.datetime.utcnow(),
                        atoms=[{'symbol': atom.symbol,
                                'position': list(atom.position),
                                'tag': atom.tag,
                                'index': atom.index,
                                'charge': atom.charge,
                                'momentum': atom.momentum.tolist(),
                                'magmom': atom.magmom} for atom in atoms],
                        pbc=atoms.pbc.tolist(),
                        info=atoms.info,
                        constraints=pickle.dumps(atoms.constraints),
                        # I would like this, but todict leaves arrays in which do not convert to json.
                        # constraints=[c.todict() for c in atoms.constraints],
                        cell=atoms.cell.tolist())

        # Convenience values
        cell = atoms.get_cell()
        if cell is not None and np.linalg.det(cell) > 0:
            d['volume'] = atoms.get_volume()

        d['mass'] = sum(atoms.get_masses())

        # fields for each symbol counts
        syms = atoms.get_chemical_symbols()
        for sym in set(syms):
            d[sym] = syms.count(sym)

        d['natoms'] = len(atoms)

        # Calculated values
        if atoms.get_calculator() is not None:
            # Need some calculator data
            calc = atoms.get_calculator()
            d['calc'] = calc.todict()

        # The rest of the data
        d.update(kwargs)

        return self.collection.insert_one(d).inserted_id

    def find(self, *args, **kwargs):
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
                                charge=atom['charge']) for atom in doc['atoms']],
                          cell=doc['cell'])

            # TODO the calculator
            yield atoms
