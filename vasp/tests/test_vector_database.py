"""Tests for the vector-enabled atomic structure database."""

import os
import tempfile
from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk, fcc111, molecule

from vasp.database import AtomMatch, MockEmbedder, VectorAtomDatabase
from vasp.database.embedders import Embedder


class TestMockEmbedder:
    """Tests for MockEmbedder."""

    def test_init_default_dim(self):
        """Test default embedding dimension."""
        embedder = MockEmbedder()
        assert embedder.dim == 128

    def test_init_custom_dim(self):
        """Test custom embedding dimension."""
        embedder = MockEmbedder(dim=64)
        assert embedder.dim == 64

    def test_embed_shape(self):
        """Test embedding output shape."""
        embedder = MockEmbedder(dim=128)
        atoms = bulk('Cu', 'fcc', a=3.6)

        embeddings = embedder.embed(atoms)

        assert embeddings.shape == (len(atoms), 128)
        assert embeddings.dtype == np.float32

    def test_embed_normalized(self):
        """Test that embeddings are unit normalized."""
        embedder = MockEmbedder(dim=128)
        atoms = bulk('Cu', 'fcc', a=3.6)

        embeddings = embedder.embed(atoms)
        norms = np.linalg.norm(embeddings, axis=1)

        np.testing.assert_allclose(norms, 1.0, rtol=1e-5)

    def test_embed_deterministic(self):
        """Test that embeddings are deterministic."""
        embedder = MockEmbedder(dim=128)
        atoms = bulk('Cu', 'fcc', a=3.6)

        emb1 = embedder.embed(atoms)
        emb2 = embedder.embed(atoms)

        np.testing.assert_array_equal(emb1, emb2)

    def test_embed_different_atoms(self):
        """Test that different atoms get different embeddings."""
        embedder = MockEmbedder(dim=128)

        cu = bulk('Cu', 'fcc', a=3.6)
        fe = bulk('Fe', 'bcc', a=2.87)

        emb_cu = embedder.embed(cu)
        emb_fe = embedder.embed(fe)

        # Embeddings should be different
        assert not np.allclose(emb_cu[0], emb_fe[0])

    def test_embed_molecule(self):
        """Test embedding a molecule."""
        embedder = MockEmbedder(dim=64)
        h2o = molecule('H2O')

        embeddings = embedder.embed(h2o)

        assert embeddings.shape == (3, 64)  # 3 atoms in H2O


class TestEmbedderInterface:
    """Test that MockEmbedder implements Embedder interface."""

    def test_is_embedder(self):
        """Test MockEmbedder is an Embedder."""
        embedder = MockEmbedder()
        assert isinstance(embedder, Embedder)

    def test_has_dim_property(self):
        """Test dim property exists."""
        embedder = MockEmbedder(dim=256)
        assert hasattr(embedder, 'dim')
        assert embedder.dim == 256

    def test_has_embed_method(self):
        """Test embed method exists."""
        embedder = MockEmbedder()
        assert hasattr(embedder, 'embed')
        assert callable(embedder.embed)


class TestVectorAtomDatabase:
    """Tests for VectorAtomDatabase."""

    @pytest.fixture
    def temp_db(self):
        """Create a temporary database."""
        with tempfile.NamedTemporaryFile(suffix='.db', delete=False) as f:
            db_path = f.name

        yield db_path

        # Cleanup
        if os.path.exists(db_path):
            os.unlink(db_path)

    @pytest.fixture
    def embedder(self):
        """Create a mock embedder."""
        return MockEmbedder(dim=64)

    @pytest.fixture
    def db(self, temp_db, embedder):
        """Create a connected database."""
        db = VectorAtomDatabase(temp_db, embedder=embedder)
        db.connect()
        yield db
        db.close()

    def test_init(self, temp_db, embedder):
        """Test database initialization."""
        db = VectorAtomDatabase(temp_db, embedder=embedder)

        assert db.db_path == Path(temp_db)
        assert db.embedder is embedder
        assert db.embedding_dim == 64

    def test_init_no_embedder(self, temp_db):
        """Test initialization without embedder."""
        db = VectorAtomDatabase(temp_db, embedding_dim=128)

        assert db.embedder is None
        assert db.embedding_dim == 128

    def test_connect(self, temp_db, embedder):
        """Test database connection."""
        db = VectorAtomDatabase(temp_db, embedder=embedder)
        db.connect()

        assert db._conn is not None
        db.close()

    def test_context_manager(self, temp_db, embedder):
        """Test context manager interface."""
        with VectorAtomDatabase(temp_db, embedder=embedder) as db:
            assert db._conn is not None

    def test_write_structure(self, db):
        """Test writing a structure."""
        atoms = bulk('Cu', 'fcc', a=3.6)
        structure_id = db.write(atoms, energy=-3.5)

        assert structure_id == 1
        assert db.count_structures() == 1
        assert db.count_atoms() == len(atoms)

    def test_write_multiple_structures(self, db):
        """Test writing multiple structures."""
        cu = bulk('Cu', 'fcc', a=3.6)
        fe = bulk('Fe', 'bcc', a=2.87)

        id1 = db.write(cu, energy=-3.5)
        id2 = db.write(fe, energy=-4.0)

        assert id1 == 1
        assert id2 == 2
        assert db.count_structures() == 2

    def test_write_with_metadata(self, db):
        """Test writing with metadata."""
        atoms = bulk('Cu', 'fcc', a=3.6)
        structure_id = db.write(
            atoms,
            energy=-3.5,
            source='DFT',
            xc='PBE',
            kpts=(8, 8, 8),
        )

        assert structure_id == 1

    def test_write_with_forces(self, db):
        """Test writing with forces."""
        atoms = bulk('Cu', 'fcc', a=3.6)
        forces = np.random.randn(len(atoms), 3)

        structure_id = db.write(atoms, forces=forces)

        assert structure_id == 1

    def test_write_without_embedder_requires_embeddings(self, temp_db):
        """Test that write fails without embedder or embeddings."""
        db = VectorAtomDatabase(temp_db)
        db.connect()

        atoms = bulk('Cu', 'fcc', a=3.6)

        with pytest.raises(ValueError, match="No embedder"):
            db.write(atoms)

        db.close()

    def test_write_with_provided_embeddings(self, temp_db):
        """Test writing with pre-computed embeddings."""
        db = VectorAtomDatabase(temp_db, embedding_dim=32)
        db.connect()

        atoms = bulk('Cu', 'fcc', a=3.6)
        embeddings = np.random.randn(len(atoms), 32).astype(np.float32)

        structure_id = db.write(atoms, embeddings=embeddings)

        assert structure_id == 1
        db.close()

    def test_get_structure(self, db):
        """Test retrieving a structure."""
        original = bulk('Cu', 'fcc', a=3.6)
        structure_id = db.write(original, energy=-3.5)

        retrieved = db.get_structure(structure_id)

        assert retrieved is not None
        assert len(retrieved) == len(original)
        assert retrieved.get_chemical_formula() == original.get_chemical_formula()
        np.testing.assert_allclose(
            retrieved.get_positions(),
            original.get_positions(),
            rtol=1e-10,
        )

    def test_get_structure_not_found(self, db):
        """Test retrieving non-existent structure."""
        result = db.get_structure(999)
        assert result is None

    def test_find_similar_atoms(self, db):
        """Test finding similar atoms."""
        # Add some structures
        cu = bulk('Cu', 'fcc', a=3.6)
        fe = bulk('Fe', 'bcc', a=2.87)

        db.write(cu, energy=-3.5)
        db.write(fe, energy=-4.0)

        # Query with embedding similar to Cu
        query_emb = db.embedder.embed(cu)[0]

        matches = db.find_similar_atoms(query_emb, k=5)

        assert len(matches) > 0
        assert isinstance(matches[0], AtomMatch)
        assert matches[0].symbol in ['Cu', 'Fe']
        # First match should be the most similar
        assert matches[0].distance <= matches[-1].distance

    def test_find_similar_atoms_filter_by_symbol(self, db):
        """Test filtering similarity search by element."""
        cu = bulk('Cu', 'fcc', a=3.6)
        fe = bulk('Fe', 'bcc', a=2.87)

        db.write(cu)
        db.write(fe)

        query_emb = db.embedder.embed(cu)[0]

        # Only search for Cu atoms
        matches = db.find_similar_atoms(query_emb, k=5, symbol='Cu')

        for match in matches:
            assert match.symbol == 'Cu'

    def test_find_similar_environments(self, db):
        """Test finding similar atomic environments."""
        cu = bulk('Cu', 'fcc', a=3.6)
        db.write(cu)

        # Find environments similar to first atom in Cu
        matches = db.find_similar_environments(cu, center_index=0, k=5)

        assert len(matches) > 0
        assert matches[0].distance < 0.1  # Should be very similar

    def test_find_similar_environments_no_embedder(self, temp_db):
        """Test that find_similar_environments fails without embedder."""
        db = VectorAtomDatabase(temp_db, embedding_dim=64)
        db.connect()

        atoms = bulk('Cu', 'fcc', a=3.6)

        with pytest.raises(ValueError, match="Embedder required"):
            db.find_similar_environments(atoms, 0)

        db.close()

    def test_count_structures(self, db):
        """Test counting structures."""
        assert db.count_structures() == 0

        db.write(bulk('Cu', 'fcc', a=3.6))
        assert db.count_structures() == 1

        db.write(bulk('Fe', 'bcc', a=2.87))
        assert db.count_structures() == 2

    def test_count_atoms(self, db):
        """Test counting atoms."""
        assert db.count_atoms() == 0

        cu = bulk('Cu', 'fcc', a=3.6)
        db.write(cu)
        assert db.count_atoms() == len(cu)

    def test_iter_structures(self, db):
        """Test iterating over structures."""
        cu = bulk('Cu', 'fcc', a=3.6)
        fe = bulk('Fe', 'bcc', a=2.87)

        db.write(cu)
        db.write(fe)

        structures = list(db.iter_structures())

        assert len(structures) == 2
        assert structures[0][0] == 1  # First ID
        assert structures[1][0] == 2  # Second ID

    def test_get_atom_embedding(self, db):
        """Test retrieving atom embeddings."""
        cu = bulk('Cu', 'fcc', a=3.6)
        db.write(cu)

        emb = db.get_atom_embedding(structure_id=1, atom_index=0)

        assert emb is not None
        assert emb.shape == (64,)
        assert emb.dtype == np.float32

    def test_get_atom_embedding_not_found(self, db):
        """Test retrieving non-existent embedding."""
        emb = db.get_atom_embedding(structure_id=999, atom_index=0)
        assert emb is None

    def test_delete_structure(self, db):
        """Test deleting a structure."""
        cu = bulk('Cu', 'fcc', a=3.6)
        structure_id = db.write(cu)

        assert db.count_structures() == 1
        assert db.count_atoms() == len(cu)

        result = db.delete_structure(structure_id)

        assert result is True
        assert db.count_structures() == 0
        assert db.count_atoms() == 0

    def test_delete_structure_not_found(self, db):
        """Test deleting non-existent structure."""
        result = db.delete_structure(999)
        assert result is False


class TestVectorDatabaseSlab:
    """Test with slab structures (more realistic use case)."""

    @pytest.fixture
    def temp_db(self):
        """Create a temporary database."""
        with tempfile.NamedTemporaryFile(suffix='.db', delete=False) as f:
            db_path = f.name

        yield db_path

        if os.path.exists(db_path):
            os.unlink(db_path)

    @pytest.fixture
    def db(self, temp_db):
        """Create database with embedder."""
        embedder = MockEmbedder(dim=128)
        db = VectorAtomDatabase(temp_db, embedder=embedder)
        db.connect()
        yield db
        db.close()

    def test_slab_structure(self, db):
        """Test storing a slab structure."""
        slab = fcc111('Pt', size=(2, 2, 4), vacuum=10.0)

        structure_id = db.write(slab, energy=-50.0)

        assert structure_id == 1
        assert db.count_atoms() == len(slab)

    def test_find_surface_atoms(self, db):
        """Test finding similar surface atomic environments."""
        # Create two similar slabs
        pt_slab = fcc111('Pt', size=(2, 2, 4), vacuum=10.0)
        pd_slab = fcc111('Pd', size=(2, 2, 4), vacuum=10.0)

        db.write(pt_slab)
        db.write(pd_slab)

        # Find atoms similar to top Pt surface atom
        top_atom_idx = len(pt_slab) - 1  # Approximate top atom

        matches = db.find_similar_environments(
            pt_slab,
            center_index=top_atom_idx,
            k=10,
        )

        assert len(matches) > 0
        # Top of slab should match other surface atoms

    def test_filter_by_element_in_slab(self, db):
        """Test filtering by element in multi-element slab."""
        from ase.build import add_adsorbate

        slab = fcc111('Pt', size=(2, 2, 4), vacuum=10.0)
        add_adsorbate(slab, 'O', height=1.5, position='fcc')

        db.write(slab)

        # Search only for O atoms
        query_emb = np.random.randn(128).astype(np.float32)
        matches = db.find_similar_atoms(query_emb, k=10, symbol='O')

        for match in matches:
            assert match.symbol == 'O'


class TestAtomMatch:
    """Tests for AtomMatch dataclass."""

    def test_creation(self):
        """Test AtomMatch creation."""
        match = AtomMatch(
            structure_id=1,
            atom_index=0,
            symbol='Cu',
            position=np.array([0.0, 0.0, 0.0]),
            distance=0.1,
            atomic_number=29,
        )

        assert match.structure_id == 1
        assert match.atom_index == 0
        assert match.symbol == 'Cu'
        assert match.distance == 0.1
        assert match.atomic_number == 29

    def test_position_array(self):
        """Test position is accessible as array."""
        match = AtomMatch(
            structure_id=1,
            atom_index=0,
            symbol='Fe',
            position=np.array([1.0, 2.0, 3.0]),
            distance=0.5,
        )

        np.testing.assert_array_equal(match.position, [1.0, 2.0, 3.0])
