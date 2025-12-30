"""Vector-enabled atomic structure database.

This module provides an ASE-compatible database with per-atom vector embeddings
stored in libSQL for similarity search.

Example:
    >>> from vasp.database import VectorAtomDatabase, MockEmbedder
    >>> from ase.build import bulk
    >>>
    >>> db = VectorAtomDatabase("atoms.db", embedder=MockEmbedder(dim=128))
    >>> await db.connect()
    >>>
    >>> atoms = bulk('Cu', 'fcc', a=3.6)
    >>> structure_id = await db.write(atoms, energy=-3.5)
    >>>
    >>> # Find similar atomic environments
    >>> similar = await db.find_similar_atoms(query_embedding, k=10)
"""

from .embedders import Embedder, MockEmbedder
from .vector_db import AtomMatch, StructureRecord, VectorAtomDatabase

__all__ = [
    'VectorAtomDatabase',
    'AtomMatch',
    'StructureRecord',
    'Embedder',
    'MockEmbedder',
]


def get_fairchem_embedder():
    """Get FairChem embedder (requires fairchem-core).

    Returns:
        FairChemEmbedder instance.

    Raises:
        ImportError: If fairchem-core is not installed.
    """
    from .embedders import FairChemEmbedder
    return FairChemEmbedder
