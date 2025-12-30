"""VASP calculation recipes for workflow integration.

Provides pre-configured calculation workflows compatible with:
- Prefect
- Dask
- Parsl
- Covalent
- Jobflow

Recipes follow the quacc pattern with @job and @flow decorators.
"""

from .core import double_relax_flow, relax_job, static_job
from .decorators import flow, job, subflow
from .phonons import phonon_flow
from .slabs import bulk_to_slabs_flow, slab_relax_job, slab_static_job

__all__ = [
    # Core jobs
    'static_job',
    'relax_job',
    'double_relax_flow',
    # Slab jobs
    'slab_static_job',
    'slab_relax_job',
    'bulk_to_slabs_flow',
    # Phonon workflows
    'phonon_flow',
    # Decorators
    'job',
    'flow',
    'subflow',
]
