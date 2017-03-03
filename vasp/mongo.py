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
                        constraints=atoms.constraints,
                        cell=atoms.cell.tolist())

        # Calculated values
        if atoms.get_calculator() is not None:
            # Need some calculator data
            calc = atoms.get_calculator()
            d['calc'] = calc.todict()

            for property in calc.implemented_properties:
                try:
                    d.update({property: getattr(calc, property)})
                except:
                    print('Failed to add {}.'.format(property))

            if 'forces' in calc.implemented_properties:
                d['fmax'] = max(np.abs(atoms.get_forces().flatten()))

            if 'stress' in calc.implemented_properties:
                d['smax'] = max(np.abs(atoms.get_stress().flatten()))

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
