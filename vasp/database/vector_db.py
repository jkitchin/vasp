"""Vector-enabled atomic structure database using libSQL.

This module provides VectorAtomDatabase, an ASE-compatible database that stores
per-atom embeddings for similarity search using libSQL's native vector support.
"""

from __future__ import annotations

import json
import sqlite3
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from ase import Atoms

    from .embedders import Embedder


@dataclass
class AtomMatch:
    """Result from atom similarity search."""

    structure_id: int
    atom_index: int
    symbol: str
    position: np.ndarray
    distance: float
    atomic_number: int = 0


@dataclass
class StructureRecord:
    """Record of a structure in the database."""

    id: int
    formula: str
    natoms: int
    energy: float | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


def _array_to_vector_string(arr: np.ndarray) -> str:
    """Convert numpy array to libSQL vector string format."""
    flat = arr.flatten().astype(np.float32)
    return "[" + ",".join(f"{x:.8g}" for x in flat) + "]"


def _vector_string_to_array(s: str) -> np.ndarray:
    """Convert libSQL vector string to numpy array."""
    # Remove brackets and split
    s = s.strip("[]")
    return np.array([float(x) for x in s.split(",")], dtype=np.float32)


class VectorAtomDatabase:
    """ASE-compatible database with per-atom vector embeddings in libSQL.

    This database stores atomic structures with per-atom embeddings generated
    by a configurable embedder (e.g., FairChem models). Embeddings are stored
    in libSQL with vector indexing for efficient similarity search.

    Features:
    - Per-atom embeddings for fine-grained similarity search
    - Find atoms with similar local environments
    - Filter by element type
    - ASE-compatible structure storage

    Example:
        >>> from vasp.database import VectorAtomDatabase, MockEmbedder
        >>> from ase.build import bulk
        >>>
        >>> # Create database with mock embedder
        >>> embedder = MockEmbedder(dim=128)
        >>> db = VectorAtomDatabase("atoms.db", embedder=embedder)
        >>> db.connect()
        >>>
        >>> # Store a structure
        >>> atoms = bulk('Cu', 'fcc', a=3.6)
        >>> structure_id = db.write(atoms, energy=-3.5)
        >>>
        >>> # Find similar atomic environments
        >>> similar = db.find_similar_atoms(query_embedding, k=10)

    Note:
        Requires libsql-experimental-python or libsql-client for full
        vector search support. Falls back to exact search without it.
    """

    def __init__(
        self,
        db_path: str | Path,
        embedder: Embedder | None = None,
        embedding_dim: int = 128,
    ):
        """Initialize vector database.

        Args:
            db_path: Path to SQLite/libSQL database file.
            embedder: Embedder instance for generating per-atom embeddings.
                     If None, embeddings must be provided manually.
            embedding_dim: Dimension of embeddings (used if embedder is None).
        """
        self.db_path = Path(db_path)
        self.embedder = embedder
        self._embedding_dim = embedding_dim if embedder is None else embedder.dim
        self._conn: sqlite3.Connection | Any | None = None
        self._has_vector_support = False

    @property
    def embedding_dim(self) -> int:
        """Get embedding dimension."""
        return self._embedding_dim

    def connect(self) -> None:
        """Connect to database and create tables."""
        # Try libSQL first for vector support
        try:
            import libsql_experimental as libsql

            self._conn = libsql.connect(str(self.db_path))
            self._has_vector_support = True
        except ImportError:
            # Fall back to standard sqlite3
            self._conn = sqlite3.connect(str(self.db_path))
            self._has_vector_support = False

        self._create_tables()

    def close(self) -> None:
        """Close database connection."""
        if self._conn:
            self._conn.close()
            self._conn = None

    def __enter__(self):
        self.connect()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _create_tables(self) -> None:
        """Create database tables and indexes."""
        cursor = self._conn.cursor()

        # Structures table
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS structures (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                formula TEXT NOT NULL,
                natoms INTEGER NOT NULL,
                cell BLOB,
                pbc BLOB,
                positions BLOB,
                numbers BLOB,
                energy REAL,
                forces BLOB,
                stress BLOB,
                metadata TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """
        )

        # Atoms table with embeddings
        if self._has_vector_support:
            cursor.execute(
                f"""
                CREATE TABLE IF NOT EXISTS atoms (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    structure_id INTEGER NOT NULL,
                    atom_index INTEGER NOT NULL,
                    symbol TEXT NOT NULL,
                    atomic_number INTEGER NOT NULL,
                    x REAL NOT NULL,
                    y REAL NOT NULL,
                    z REAL NOT NULL,
                    fx REAL,
                    fy REAL,
                    fz REAL,
                    embedding F32_BLOB({self._embedding_dim}),
                    FOREIGN KEY (structure_id) REFERENCES structures(id),
                    UNIQUE(structure_id, atom_index)
                )
            """
            )

            # Create vector index
            try:
                cursor.execute(
                    """
                    CREATE INDEX IF NOT EXISTS atoms_emb_idx ON atoms (
                        libsql_vector_idx(embedding, 'metric=cosine')
                    )
                """
                )
            except Exception:
                # Vector index creation might fail on some libSQL versions
                pass
        else:
            # Standard SQLite without vector support
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS atoms (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    structure_id INTEGER NOT NULL,
                    atom_index INTEGER NOT NULL,
                    symbol TEXT NOT NULL,
                    atomic_number INTEGER NOT NULL,
                    x REAL NOT NULL,
                    y REAL NOT NULL,
                    z REAL NOT NULL,
                    fx REAL,
                    fy REAL,
                    fz REAL,
                    embedding BLOB,
                    FOREIGN KEY (structure_id) REFERENCES structures(id),
                    UNIQUE(structure_id, atom_index)
                )
            """
            )

        # Index for element queries
        cursor.execute("CREATE INDEX IF NOT EXISTS atoms_symbol_idx ON atoms(symbol)")
        cursor.execute("CREATE INDEX IF NOT EXISTS atoms_structure_idx ON atoms(structure_id)")

        self._conn.commit()

    def write(
        self,
        atoms: Atoms,
        energy: float | None = None,
        forces: np.ndarray | None = None,
        stress: np.ndarray | None = None,
        embeddings: np.ndarray | None = None,
        **metadata,
    ) -> int:
        """Write structure with per-atom embeddings to database.

        Args:
            atoms: ASE Atoms object to store.
            energy: Potential energy in eV.
            forces: Forces array of shape (natoms, 3) in eV/Angstrom.
            stress: Stress tensor (6,) in Voigt notation.
            embeddings: Pre-computed embeddings of shape (natoms, dim).
                       If None, uses self.embedder to generate them.
            **metadata: Additional key-value pairs to store.

        Returns:
            Structure ID in database.

        Raises:
            ValueError: If no embedder and no embeddings provided.
        """
        # Get embeddings
        if embeddings is None:
            if self.embedder is None:
                raise ValueError(
                    "No embedder configured and no embeddings provided. "
                    "Either pass embedder to __init__ or provide embeddings."
                )
            embeddings = self.embedder.embed(atoms)

        if embeddings.shape[0] != len(atoms):
            raise ValueError(
                f"Embeddings shape {embeddings.shape} doesn't match "
                f"number of atoms {len(atoms)}"
            )

        cursor = self._conn.cursor()

        # Insert structure
        cursor.execute(
            """INSERT INTO structures
               (formula, natoms, cell, pbc, positions, numbers, energy, forces, stress, metadata)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                atoms.get_chemical_formula(),
                len(atoms),
                atoms.cell.array.tobytes() if atoms.cell is not None else None,
                np.array(atoms.pbc, dtype=np.uint8).tobytes(),
                atoms.positions.tobytes(),
                atoms.numbers.tobytes(),
                energy,
                forces.tobytes() if forces is not None else None,
                stress.tobytes() if stress is not None else None,
                json.dumps(metadata) if metadata else None,
            ),
        )
        structure_id = cursor.lastrowid
        assert structure_id is not None

        # Insert atoms with embeddings
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        numbers = atoms.get_atomic_numbers()

        for i in range(len(atoms)):
            emb = embeddings[i]
            fx, fy, fz = forces[i] if forces is not None else (None, None, None)

            if self._has_vector_support:
                cursor.execute(
                    """INSERT INTO atoms
                       (structure_id, atom_index, symbol, atomic_number, x, y, z, fx, fy, fz, embedding)
                       VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, vector(?))""",
                    (
                        structure_id,
                        i,
                        symbols[i],
                        int(numbers[i]),
                        float(positions[i, 0]),
                        float(positions[i, 1]),
                        float(positions[i, 2]),
                        fx,
                        fy,
                        fz,
                        _array_to_vector_string(emb),
                    ),
                )
            else:
                cursor.execute(
                    """INSERT INTO atoms
                       (structure_id, atom_index, symbol, atomic_number, x, y, z, fx, fy, fz, embedding)
                       VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                    (
                        structure_id,
                        i,
                        symbols[i],
                        int(numbers[i]),
                        float(positions[i, 0]),
                        float(positions[i, 1]),
                        float(positions[i, 2]),
                        fx,
                        fy,
                        fz,
                        emb.tobytes(),
                    ),
                )

        self._conn.commit()
        return structure_id

    def get_structure(self, structure_id: int) -> Atoms | None:
        """Retrieve structure by ID.

        Args:
            structure_id: Structure ID.

        Returns:
            ASE Atoms object or None if not found.
        """
        from ase import Atoms

        cursor = self._conn.cursor()
        cursor.execute(
            "SELECT formula, natoms, cell, pbc, positions, numbers FROM structures WHERE id = ?",
            (structure_id,),
        )
        row = cursor.fetchone()

        if row is None:
            return None

        formula, natoms, cell_bytes, pbc_bytes, pos_bytes, num_bytes = row

        positions = np.frombuffer(pos_bytes, dtype=np.float64).reshape(-1, 3)
        numbers = np.frombuffer(num_bytes, dtype=np.int64)
        cell = np.frombuffer(cell_bytes, dtype=np.float64).reshape(3, 3) if cell_bytes else None
        pbc = np.frombuffer(pbc_bytes, dtype=np.uint8).astype(bool) if pbc_bytes else False

        atoms = Atoms(numbers=numbers, positions=positions, cell=cell, pbc=pbc)
        return atoms

    def find_similar_atoms(
        self,
        query_embedding: np.ndarray,
        k: int = 10,
        symbol: str | None = None,
    ) -> list[AtomMatch]:
        """Find k most similar atoms by embedding.

        Args:
            query_embedding: Query embedding vector.
            k: Number of results to return.
            symbol: Filter by element symbol (e.g., 'Fe').

        Returns:
            List of AtomMatch objects sorted by distance.
        """
        cursor = self._conn.cursor()
        query_vec = _array_to_vector_string(query_embedding)

        if self._has_vector_support:
            # Use libSQL vector search
            if symbol:
                cursor.execute(
                    """SELECT a.structure_id, a.atom_index, a.symbol,
                              a.atomic_number, a.x, a.y, a.z, v.distance
                       FROM vector_top_k('atoms_emb_idx', vector(?), ?) AS v
                       JOIN atoms a ON a.rowid = v.id
                       WHERE a.symbol = ?
                       LIMIT ?""",
                    (query_vec, k * 3, symbol, k),
                )
            else:
                cursor.execute(
                    """SELECT a.structure_id, a.atom_index, a.symbol,
                              a.atomic_number, a.x, a.y, a.z, v.distance
                       FROM vector_top_k('atoms_emb_idx', vector(?), ?) AS v
                       JOIN atoms a ON a.rowid = v.id""",
                    (query_vec, k),
                )
        else:
            # Fallback: exact search (slow for large databases)
            if symbol:
                cursor.execute(
                    "SELECT id, structure_id, atom_index, symbol, atomic_number, x, y, z, embedding FROM atoms WHERE symbol = ?",
                    (symbol,),
                )
            else:
                cursor.execute(
                    "SELECT id, structure_id, atom_index, symbol, atomic_number, x, y, z, embedding FROM atoms"
                )

            # Compute distances manually
            results = []
            query_norm = query_embedding / (np.linalg.norm(query_embedding) + 1e-8)

            for row in cursor.fetchall():
                emb_bytes = row[8]
                if emb_bytes:
                    emb = np.frombuffer(emb_bytes, dtype=np.float32)
                    emb_norm = emb / (np.linalg.norm(emb) + 1e-8)
                    # Cosine distance = 1 - cosine similarity
                    distance = 1.0 - float(np.dot(query_norm, emb_norm))
                    results.append((row, distance))

            # Sort by distance and take top k
            results.sort(key=lambda x: x[1])
            results = results[:k]

            return [
                AtomMatch(
                    structure_id=row[1],
                    atom_index=row[2],
                    symbol=row[3],
                    atomic_number=row[4],
                    position=np.array([row[5], row[6], row[7]]),
                    distance=dist,
                )
                for row, dist in results
            ]

        # Process libSQL vector search results
        matches = []
        for row in cursor.fetchall():
            matches.append(
                AtomMatch(
                    structure_id=row[0],
                    atom_index=row[1],
                    symbol=row[2],
                    atomic_number=row[3],
                    position=np.array([row[4], row[5], row[6]]),
                    distance=row[7],
                )
            )

        return matches

    def find_similar_environments(
        self,
        atoms: Atoms,
        center_index: int,
        k: int = 10,
        symbol: str | None = None,
    ) -> list[AtomMatch]:
        """Find atoms with similar local environments.

        Args:
            atoms: ASE Atoms object containing the query atom.
            center_index: Index of the query atom in atoms.
            k: Number of results to return.
            symbol: Filter by element symbol.

        Returns:
            List of AtomMatch objects.

        Raises:
            ValueError: If no embedder configured.
        """
        if self.embedder is None:
            raise ValueError("Embedder required for environment search")

        embeddings = self.embedder.embed(atoms)
        query_emb = embeddings[center_index]

        return self.find_similar_atoms(query_emb, k=k, symbol=symbol)

    def count_structures(self) -> int:
        """Count total structures in database."""
        cursor = self._conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM structures")
        return cursor.fetchone()[0]

    def count_atoms(self) -> int:
        """Count total atoms in database."""
        cursor = self._conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM atoms")
        return cursor.fetchone()[0]

    def iter_structures(self) -> Iterator[tuple[int, Atoms]]:
        """Iterate over all structures in database.

        Yields:
            Tuples of (structure_id, atoms).
        """
        cursor = self._conn.cursor()
        cursor.execute("SELECT id FROM structures ORDER BY id")

        for (structure_id,) in cursor.fetchall():
            atoms = self.get_structure(structure_id)
            if atoms is not None:
                yield structure_id, atoms

    def get_atom_embedding(self, structure_id: int, atom_index: int) -> np.ndarray | None:
        """Get embedding for a specific atom.

        Args:
            structure_id: Structure ID.
            atom_index: Atom index within structure.

        Returns:
            Embedding array or None if not found.
        """
        cursor = self._conn.cursor()
        cursor.execute(
            "SELECT embedding FROM atoms WHERE structure_id = ? AND atom_index = ?",
            (structure_id, atom_index),
        )
        row = cursor.fetchone()

        if row is None or row[0] is None:
            return None

        if self._has_vector_support:
            # Vector format - need to parse
            return _vector_string_to_array(row[0])
        else:
            return np.frombuffer(row[0], dtype=np.float32)

    def delete_structure(self, structure_id: int) -> bool:
        """Delete a structure and its atoms.

        Args:
            structure_id: Structure ID to delete.

        Returns:
            True if structure was deleted, False if not found.
        """
        cursor = self._conn.cursor()
        cursor.execute("DELETE FROM atoms WHERE structure_id = ?", (structure_id,))
        cursor.execute("DELETE FROM structures WHERE id = ?", (structure_id,))
        self._conn.commit()
        return cursor.rowcount > 0
