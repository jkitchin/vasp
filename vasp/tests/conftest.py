"""Pytest configuration and shared fixtures."""

import os
import shutil
import tempfile

import numpy as np
import pytest

from vasp import MockRunner, Vasp
from vasp.runners import JobState, MockResults


@pytest.fixture
def temp_dir():
    """Create a temporary directory for tests."""
    d = tempfile.mkdtemp()
    yield d
    shutil.rmtree(d)


@pytest.fixture
def calc_dir(temp_dir):
    """Create a calculation directory with mock input files."""
    from vasp.tests.fixtures import MOCK_INCAR, MOCK_KPOINTS, MOCK_POSCAR, MOCK_POTCAR_HEADER

    # Write mock input files
    with open(os.path.join(temp_dir, 'INCAR'), 'w') as f:
        f.write(MOCK_INCAR)

    with open(os.path.join(temp_dir, 'POSCAR'), 'w') as f:
        f.write(MOCK_POSCAR)

    with open(os.path.join(temp_dir, 'KPOINTS'), 'w') as f:
        f.write(MOCK_KPOINTS)

    with open(os.path.join(temp_dir, 'POTCAR'), 'w') as f:
        f.write(MOCK_POTCAR_HEADER)

    return temp_dir


@pytest.fixture
def complete_calc_dir(calc_dir):
    """Create a calculation directory with completed outputs."""
    from vasp.tests.fixtures import MOCK_OUTCAR_COMPLETE, MOCK_VASPRUN_XML

    with open(os.path.join(calc_dir, 'OUTCAR'), 'w') as f:
        f.write(MOCK_OUTCAR_COMPLETE)

    with open(os.path.join(calc_dir, 'vasprun.xml'), 'w') as f:
        f.write(MOCK_VASPRUN_XML)

    return calc_dir


@pytest.fixture
def si_atoms():
    """Create a silicon diamond structure."""
    from ase import Atoms

    a = 5.43
    atoms = Atoms(
        'Si2',
        scaled_positions=[(0, 0, 0), (0.25, 0.25, 0.25)],
        cell=[[a, 0, 0], [0, a, 0], [0, 0, a]],
        pbc=True
    )
    return atoms


@pytest.fixture
def mock_results():
    """Standard mock calculation results."""
    return MockResults(
        energy=-10.84274516,
        forces=np.array([
            [0.012345, -0.023456, 0.001234],
            [-0.012345, 0.023456, -0.001234]
        ]),
        stress=np.array([1.234, 1.234, 1.234, -0.123, 0.012, 0.012]),
        magmom=0.0,
        magmoms=np.array([0.0, 0.0])
    )


@pytest.fixture
def mock_runner(mock_results):
    """Create a mock runner with standard results."""
    return MockRunner(results=mock_results)


@pytest.fixture
def vasp_calc(temp_dir, si_atoms, mock_runner):
    """Create a Vasp calculator with mock runner that has completed calculation."""
    calc = Vasp(
        label=temp_dir,
        atoms=si_atoms,
        runner=mock_runner,
        xc='PBE',
        encut=400,
        kpts=(8, 8, 8)
    )
    # Run the calculation so results are available
    calc.calculate()
    return calc


@pytest.fixture
def vasp_calc_fresh(temp_dir, si_atoms, mock_runner):
    """Create a Vasp calculator with mock runner (not yet calculated)."""
    calc = Vasp(
        label=temp_dir,
        atoms=si_atoms,
        runner=mock_runner,
        xc='PBE',
        encut=400,
        kpts=(8, 8, 8)
    )
    return calc


@pytest.fixture
def async_mock_runner(mock_results):
    """Create a mock runner that simulates async workflow."""
    return MockRunner(
        results=mock_results,
        state_sequence=[JobState.SUBMITTED, JobState.QUEUED, JobState.RUNNING, JobState.COMPLETE]
    )
