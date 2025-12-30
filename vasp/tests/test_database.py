"""Tests for ASE database integration with Vasp calculator.

This module tests that Vasp calculators work correctly with ASE's
database functionality, ensuring atoms and calculation data can be
stored and retrieved.
"""

import os

import numpy as np
import pytest
from ase import Atoms
from ase.db import connect

from vasp import Vasp
from vasp.runners import MockResults, MockRunner


class TestDatabaseWrite:
    """Test writing atoms and calculator data to ASE database."""

    def test_write_atoms_to_database(self, temp_dir, si_atoms, mock_runner):
        """Test that atoms with calculator can be written to database."""
        db_path = os.path.join(temp_dir, 'test.db')

        calc = Vasp(
            label=temp_dir,
            atoms=si_atoms,
            runner=mock_runner,
            xc='PBE',
            encut=400,
            kpts=(8, 8, 8),
        )
        calc.calculate()

        # Write atoms to database
        with connect(db_path) as db:
            db.write(si_atoms, calculator_name='vasp')

        # Verify database entry exists
        with connect(db_path) as db:
            assert db.count() == 1

    def test_write_with_key_value_pairs(self, temp_dir, si_atoms, mock_runner):
        """Test writing atoms with custom key-value pairs."""
        db_path = os.path.join(temp_dir, 'test.db')

        calc = Vasp(
            label=temp_dir,
            atoms=si_atoms,
            runner=mock_runner,
            sigma=0.01,
        )
        calc.calculate()

        # Write with key-value pairs (similar to legacy test)
        with connect(db_path) as db:
            db.write(
                si_atoms,
                key_value_pairs={'xc': 'PBE', 'project': 'test'},
                data={'parameters': calc.parameters},
            )

        # Verify data can be retrieved
        with connect(db_path) as db:
            row = db.get(id=1)
            assert row.key_value_pairs['xc'] == 'PBE'
            assert row.key_value_pairs['project'] == 'test'
            assert row.data['parameters']['sigma'] == 0.01

    def test_write_calculation_results(self, temp_dir, si_atoms, mock_runner):
        """Test storing calculation results in database."""
        db_path = os.path.join(temp_dir, 'test.db')

        calc = Vasp(
            label=temp_dir,
            atoms=si_atoms,
            runner=mock_runner,
        )
        calc.calculate()

        energy = calc.results['energy']
        forces = calc.results['forces']

        # Store results
        with connect(db_path) as db:
            db.write(
                si_atoms,
                data={
                    'energy': energy,
                    'forces': forces.tolist(),
                    'path': calc.directory,
                },
            )

        # Verify results retrieval
        with connect(db_path) as db:
            row = db.get(id=1)
            assert row.data['energy'] == pytest.approx(-10.84274516)
            assert len(row.data['forces']) == 2

    def test_write_magnetic_atoms(self, temp_dir):
        """Test writing atoms with magnetic moments (like legacy CO2 test)."""
        db_path = os.path.join(temp_dir, 'test.db')

        # Create CO2 molecule with magnetic moments (similar to legacy test)
        atoms = Atoms(
            symbols=['O', 'C', 'O'],
            positions=[
                [4, 5, 5],
                [5, 5, 5],
                [6, 5, 5],
            ],
            cell=(10, 10, 10),
            pbc=True,
        )
        atoms.set_initial_magnetic_moments([1.0, 2.0, 3.0])

        mock = MockResults(
            energy=-15.0,
            forces=np.zeros((3, 3)),
            magmom=6.0,
            magmoms=np.array([1.0, 2.0, 3.0]),
        )
        runner = MockRunner(results=mock)

        calc = Vasp(
            label=temp_dir,
            atoms=atoms,
            runner=runner,
            sigma=0.01,
            ispin=2,
        )
        calc.calculate()

        # Write to database with parameters
        with connect(db_path) as db:
            db.write(
                atoms,
                data={'parameters': calc.parameters},
            )

        # Verify parameters stored correctly (key assertion from legacy test)
        with connect(db_path) as db:
            data = db.get(id=1).data
            assert data['parameters']['sigma'] == 0.01
            assert data['parameters']['ispin'] == 2


class TestDatabaseRead:
    """Test reading atoms and data from ASE database."""

    def test_read_atoms_from_database(self, temp_dir, si_atoms, mock_runner):
        """Test reading atoms back from database."""
        db_path = os.path.join(temp_dir, 'test.db')

        calc = Vasp(
            label=temp_dir,
            atoms=si_atoms,
            runner=mock_runner,
        )
        calc.calculate()

        # Write atoms
        with connect(db_path) as db:
            db.write(si_atoms)

        # Read atoms back
        with connect(db_path) as db:
            read_atoms = db.get_atoms(id=1)

        # Verify atoms match
        assert len(read_atoms) == len(si_atoms)
        assert read_atoms.get_chemical_symbols() == si_atoms.get_chemical_symbols()
        np.testing.assert_allclose(
            read_atoms.get_positions(),
            si_atoms.get_positions(),
            rtol=1e-10,
        )

    def test_read_multiple_entries(self, temp_dir, si_atoms, mock_runner):
        """Test reading multiple entries from database."""
        db_path = os.path.join(temp_dir, 'test.db')

        calc = Vasp(
            label=temp_dir,
            atoms=si_atoms,
            runner=mock_runner,
        )
        calc.calculate()

        # Write multiple entries with different energies
        energies = [-10.0, -11.0, -12.0]
        with connect(db_path) as db:
            for e in energies:
                db.write(si_atoms, data={'energy': e})

        # Read all entries
        with connect(db_path) as db:
            assert db.count() == 3
            for i, row in enumerate(db.select()):
                assert row.data['energy'] == energies[i]


class TestDatabaseQuery:
    """Test querying ASE database with calculator data."""

    def test_query_by_formula(self, temp_dir, si_atoms, mock_runner):
        """Test querying database by chemical formula."""
        db_path = os.path.join(temp_dir, 'test.db')

        calc = Vasp(
            label=temp_dir,
            atoms=si_atoms,
            runner=mock_runner,
        )
        calc.calculate()

        with connect(db_path) as db:
            db.write(si_atoms, key_value_pairs={'xc': 'PBE'})

        # Query by formula
        with connect(db_path) as db:
            results = list(db.select(formula='Si2'))
            assert len(results) == 1
            assert results[0].formula == 'Si2'

    def test_query_by_key_value(self, temp_dir, si_atoms, mock_runner):
        """Test querying by key-value pairs."""
        db_path = os.path.join(temp_dir, 'test.db')

        calc = Vasp(
            label=temp_dir,
            atoms=si_atoms,
            runner=mock_runner,
        )
        calc.calculate()

        # Write entries with different XC functionals
        with connect(db_path) as db:
            db.write(si_atoms, key_value_pairs={'xc': 'PBE', 'encut': 400})
            db.write(si_atoms, key_value_pairs={'xc': 'LDA', 'encut': 400})
            db.write(si_atoms, key_value_pairs={'xc': 'PBE', 'encut': 500})

        # Query by XC
        with connect(db_path) as db:
            pbe_results = list(db.select(xc='PBE'))
            assert len(pbe_results) == 2

            # Query by multiple criteria
            high_encut = list(db.select(xc='PBE', encut=500))
            assert len(high_encut) == 1


class TestDatabaseAppend:
    """Test appending to existing database entries."""

    def test_append_calculation(self, temp_dir, si_atoms, mock_runner):
        """Test appending new calculations to database."""
        db_path = os.path.join(temp_dir, 'test.db')

        calc = Vasp(
            label=temp_dir,
            atoms=si_atoms,
            runner=mock_runner,
        )
        calc.calculate()

        # First write
        with connect(db_path) as db:
            db.write(si_atoms, data={'step': 1})

        # Append more data
        with connect(db_path) as db:
            db.write(si_atoms, data={'step': 2})

        with connect(db_path) as db:
            assert db.count() == 2

    def test_overwrite_entry(self, temp_dir, si_atoms, mock_runner):
        """Test overwriting database entry (like legacy append=False)."""
        db_path = os.path.join(temp_dir, 'test.db')

        calc = Vasp(
            label=temp_dir,
            atoms=si_atoms,
            runner=mock_runner,
        )
        calc.calculate()

        # Write initial entry
        with connect(db_path) as db:
            row_id = db.write(si_atoms, data={'version': 1})

        # Overwrite by deleting and rewriting
        with connect(db_path) as db:
            db.delete([row_id])
            db.write(si_atoms, data={'version': 2})

        with connect(db_path) as db:
            assert db.count() == 1
            row = db.get(id=2)  # New ID after delete
            assert row.data['version'] == 2
