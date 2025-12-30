"""Workflow decorators for different execution backends.

Provides unified @job, @flow, and @subflow decorators that work with:
- No workflow engine (direct execution)
- Prefect
- Dask
- Parsl
- Covalent

The workflow engine is selected via the VASP_WORKFLOW_ENGINE environment
variable or can be set programmatically.
"""

from __future__ import annotations

import functools
import os
from typing import Any, Callable, TypeVar

F = TypeVar('F', bound=Callable[..., Any])

# Detect available workflow engines
WORKFLOW_ENGINE = os.environ.get('VASP_WORKFLOW_ENGINE', 'none').lower()

_prefect_available = False
_dask_available = False
_parsl_available = False
_covalent_available = False

try:
    import prefect  # noqa: F401
    _prefect_available = True
except ImportError:
    pass

try:
    import dask  # noqa: F401
    _dask_available = True
except ImportError:
    pass

try:
    import parsl  # noqa: F401
    _parsl_available = True
except ImportError:
    pass

try:
    import covalent  # noqa: F401
    _covalent_available = True
except ImportError:
    pass


def _identity_decorator(func: F) -> F:
    """No-op decorator for direct execution."""
    return func


def job(
    func: F | None = None,
    *,
    name: str | None = None,
    retries: int = 0,
    retry_delay_seconds: float = 0,
    **kwargs
) -> F | Callable[[F], F]:
    """Decorate a function as a workflow job.

    This decorator adapts to the configured workflow engine:
    - none: Direct function execution
    - prefect: Wraps with @prefect.task
    - dask: Wraps with @dask.delayed
    - parsl: Wraps with @parsl.python_app
    - covalent: Wraps with @covalent.electron

    Args:
        func: Function to decorate.
        name: Job name (for workflow UIs).
        retries: Number of retry attempts.
        retry_delay_seconds: Delay between retries.
        **kwargs: Additional backend-specific arguments.

    Returns:
        Decorated function.

    Example:
        >>> @job
        ... def my_calculation(atoms):
        ...     calc = Vasp(...)
        ...     return calc.get_potential_energy()
    """
    def decorator(fn: F) -> F:
        if WORKFLOW_ENGINE == 'prefect' and _prefect_available:
            from prefect import task
            return task(
                fn,
                name=name or fn.__name__,
                retries=retries,
                retry_delay_seconds=retry_delay_seconds,
                **kwargs
            )

        elif WORKFLOW_ENGINE == 'dask' and _dask_available:
            from dask import delayed
            @functools.wraps(fn)
            def wrapper(*args, **kw):
                return delayed(fn)(*args, **kw)
            return wrapper

        elif WORKFLOW_ENGINE == 'parsl' and _parsl_available:
            from parsl import python_app
            return python_app(fn)

        elif WORKFLOW_ENGINE == 'covalent' and _covalent_available:
            import covalent as ct
            return ct.electron(fn)

        else:
            # Direct execution - no wrapping
            return fn

    if func is not None:
        return decorator(func)
    return decorator


def flow(
    func: F | None = None,
    *,
    name: str | None = None,
    **kwargs
) -> F | Callable[[F], F]:
    """Decorate a function as a workflow flow.

    A flow orchestrates multiple jobs together.

    Args:
        func: Function to decorate.
        name: Flow name.
        **kwargs: Backend-specific arguments.

    Returns:
        Decorated function.

    Example:
        >>> @flow
        ... def relax_then_static(atoms):
        ...     relaxed = relax_job(atoms)
        ...     return static_job(relaxed.atoms)
    """
    def decorator(fn: F) -> F:
        if WORKFLOW_ENGINE == 'prefect' and _prefect_available:
            from prefect import flow as prefect_flow
            return prefect_flow(fn, name=name or fn.__name__, **kwargs)

        elif WORKFLOW_ENGINE == 'covalent' and _covalent_available:
            import covalent as ct
            return ct.lattice(fn)

        else:
            # Direct execution or other engines
            return fn

    if func is not None:
        return decorator(func)
    return decorator


def subflow(
    func: F | None = None,
    *,
    name: str | None = None,
    **kwargs
) -> F | Callable[[F], F]:
    """Decorate a function as a subflow (nested flow).

    Args:
        func: Function to decorate.
        name: Subflow name.
        **kwargs: Backend-specific arguments.

    Returns:
        Decorated function.
    """
    # For most engines, subflow == flow
    return flow(func, name=name, **kwargs)


def set_workflow_engine(engine: str) -> None:
    """Set the workflow engine programmatically.

    Args:
        engine: Engine name ('none', 'prefect', 'dask', 'parsl', 'covalent').
    """
    global WORKFLOW_ENGINE
    valid = {'none', 'prefect', 'dask', 'parsl', 'covalent'}
    if engine.lower() not in valid:
        raise ValueError(f"Unknown engine '{engine}'. Valid: {valid}")
    WORKFLOW_ENGINE = engine.lower()


def get_workflow_engine() -> str:
    """Get the current workflow engine name."""
    return WORKFLOW_ENGINE
