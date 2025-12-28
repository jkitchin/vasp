"""Embedder interfaces for generating per-atom vector representations.

This module provides:
- Abstract Embedder base class
- MockEmbedder for testing
- FairChemEmbedder for production use with FairChem models
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ase import Atoms


class Embedder(ABC):
    """Abstract base class for atom embedders."""

    @property
    @abstractmethod
    def dim(self) -> int:
        """Embedding dimension."""
        ...

    @abstractmethod
    def embed(self, atoms: Atoms) -> np.ndarray:
        """Generate per-atom embeddings.

        Args:
            atoms: ASE Atoms object.

        Returns:
            Array of shape (natoms, dim) with per-atom embeddings.
        """
        ...


class MockEmbedder(Embedder):
    """Mock embedder for testing.

    Generates deterministic embeddings based on atomic properties.
    """

    def __init__(self, dim: int = 128):
        """Initialize mock embedder.

        Args:
            dim: Embedding dimension.
        """
        self._dim = dim

    @property
    def dim(self) -> int:
        return self._dim

    def embed(self, atoms: Atoms) -> np.ndarray:
        """Generate mock per-atom embeddings.

        Creates deterministic embeddings based on:
        - Atomic number
        - Position (normalized)
        - Random but reproducible component based on index

        Args:
            atoms: ASE Atoms object.

        Returns:
            Array of shape (natoms, dim).
        """
        n_atoms = len(atoms)
        embeddings = np.zeros((n_atoms, self._dim), dtype=np.float32)

        positions = atoms.get_positions()
        numbers = atoms.get_atomic_numbers()

        # Normalize positions if cell exists
        if atoms.cell.volume > 0:
            pos_norm = positions / np.cbrt(atoms.cell.volume)
        else:
            pos_norm = positions / (np.max(np.abs(positions)) + 1e-8)

        for i in range(n_atoms):
            # Seed based on atom index and atomic number for reproducibility
            rng = np.random.RandomState(seed=numbers[i] * 1000 + i)

            # First few dimensions: atomic properties
            embeddings[i, 0] = numbers[i] / 100.0  # Normalized atomic number
            embeddings[i, 1:4] = pos_norm[i]  # Normalized position

            # Rest: deterministic pseudo-random based on element
            embeddings[i, 4:] = rng.randn(self._dim - 4) * 0.1

        # Normalize to unit vectors
        norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
        embeddings = embeddings / (norms + 1e-8)

        return embeddings


class FairChemEmbedder(Embedder):
    """Embedder using FairChem pretrained models.

    Uses the backbone of FairChem models (e.g., UMA) to generate
    per-atom embeddings that capture local chemical environments.

    Requires: pip install fairchem-core torch

    Example:
        >>> embedder = FairChemEmbedder(model_name="uma-s-1")
        >>> embeddings = embedder.embed(atoms)  # Shape: (natoms, 256)
    """

    def __init__(
        self,
        model_name: str = "uma-s-1",
        device: str | None = None,
    ):
        """Initialize FairChem embedder.

        Args:
            model_name: Name of pretrained model to use.
            device: Device to run model on ('cuda', 'cpu', or None for auto).
        """
        try:
            import torch
            from fairchem.core.models import load_model
        except ImportError as e:
            raise ImportError(
                "FairChemEmbedder requires fairchem-core and torch. "
                "Install with: pip install fairchem-core torch"
            ) from e

        self.model_name = model_name

        # Set device
        if device is None:
            self._device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        else:
            self._device = torch.device(device)

        # Load model and extract backbone
        self.model = load_model(model_name)
        self.backbone = self.model.backbone
        self.backbone.to(self._device)
        self.backbone.eval()

        # Determine embedding dimension from model config
        self._dim = self._get_embedding_dim()

    def _get_embedding_dim(self) -> int:
        """Determine embedding dimension from model."""
        # Try to get from model config, fallback to common default
        if hasattr(self.backbone, 'hidden_channels'):
            return self.backbone.hidden_channels
        if hasattr(self.backbone, 'emb_size'):
            return self.backbone.emb_size
        # Default for UMA models
        return 256

    @property
    def dim(self) -> int:
        return self._dim

    def embed(self, atoms: Atoms) -> np.ndarray:
        """Generate per-atom embeddings using FairChem backbone.

        Args:
            atoms: ASE Atoms object.

        Returns:
            Array of shape (natoms, dim) with per-atom embeddings.
        """
        import torch
        from fairchem.core.datasets import atoms_to_data

        # Convert ASE Atoms to FairChem AtomicData
        data = atoms_to_data(atoms)
        data = data.to(self._device)

        with torch.no_grad():
            # Forward through backbone to get embeddings
            emb_dict = self.backbone(data)

            # Extract node/atom embeddings
            # Different models may use different keys
            for key in ['node_embedding', 'h', 'x', 'node_features']:
                if key in emb_dict:
                    atom_emb = emb_dict[key]
                    break
            else:
                raise KeyError(
                    f"Could not find atom embeddings in backbone output. "
                    f"Available keys: {list(emb_dict.keys())}"
                )

        return atom_emb.cpu().numpy().astype(np.float32)
