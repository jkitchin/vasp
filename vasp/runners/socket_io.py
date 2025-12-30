"""Socket I/O for remote VASP communication.

This module implements the i-PI protocol for socket-based communication
with VASP or other DFT codes. This allows separation of the driver
(optimization algorithm) from the compute backend (VASP).

The i-PI protocol is a simple binary protocol:
1. Server sends STATUS request
2. Client responds with current state (READY, NEEDINIT, HAVEDATA)
3. Server sends POSDATA with cell and positions
4. Client runs calculation and responds with FORCEREADY
5. Server sends GETFORCE request
6. Client sends energy, forces, and virial

References:
    i-PI: https://ipi-code.org/
    ASE SocketIOCalculator: ase.calculators.socketio

Inspired by:
    VaspInteractive by the Ulissi Group at CMU
    https://github.com/ulissigroup/vasp-interactive
"""

from __future__ import annotations

import socket
import struct
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ase import Atoms

    from vasp.runners.interactive import InteractiveRunner


# =============================================================================
# Protocol Constants
# =============================================================================

# Message types (12 bytes, padded with spaces)
MSG_STATUS = b"STATUS      "
MSG_READY = b"READY       "
MSG_NEEDINIT = b"NEEDINIT    "
MSG_HAVEDATA = b"HAVEDATA    "
MSG_INIT = b"INIT        "
MSG_POSDATA = b"POSDATA     "
MSG_GETFORCE = b"GETFORCE    "
MSG_FORCEREADY = b"FORCEREADY  "
MSG_EXIT = b"EXIT        "

# Header size
HEADER_SIZE = 12

# Bohr to Angstrom conversion
BOHR_TO_ANGSTROM = 0.52917721067
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM

# Hartree to eV conversion
HARTREE_TO_EV = 27.211386245988
EV_TO_HARTREE = 1.0 / HARTREE_TO_EV


class SocketState(Enum):
    """State of the socket client."""

    DISCONNECTED = "disconnected"
    CONNECTED = "connected"
    READY = "ready"
    BUSY = "busy"


@dataclass
class SocketConfig:
    """Configuration for socket connection."""

    host: str = "localhost"
    port: int = 31415
    unix_socket: str | None = None
    timeout: float = 3600.0


# =============================================================================
# Socket Server (Driver Side)
# =============================================================================


class SocketServer:
    """Socket server for driving VASP calculations.

    This server implements the driver side of the i-PI protocol.
    It sends atomic positions to a VASP client and receives
    energies and forces.

    Args:
        config: Socket configuration.

    Example:
        >>> server = SocketServer(SocketConfig(port=31415))
        >>> server.start()
        >>> server.wait_for_client()
        >>>
        >>> # Send positions, get forces
        >>> results = server.calculate(atoms)
        >>>
        >>> server.close()
    """

    def __init__(self, config: SocketConfig | None = None):
        self.config = config or SocketConfig()
        self._socket: socket.socket | None = None
        self._client: socket.socket | None = None
        self._initialized = False

    def start(self) -> None:
        """Start the socket server."""
        if self.config.unix_socket:
            self._socket = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            self._socket.bind(self.config.unix_socket)
        else:
            self._socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self._socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            self._socket.bind((self.config.host, self.config.port))

        self._socket.listen(1)
        self._socket.settimeout(self.config.timeout)

    def wait_for_client(self) -> None:
        """Wait for a client to connect."""
        if self._socket is None:
            raise RuntimeError("Server not started")

        self._client, addr = self._socket.accept()
        self._client.settimeout(self.config.timeout)
        self._initialized = False

    def calculate(self, atoms: Atoms) -> dict:
        """Send positions and receive energy/forces.

        Args:
            atoms: Atomic structure with positions.

        Returns:
            Dict with 'energy', 'forces', and optionally 'stress'.
        """
        if self._client is None:
            raise RuntimeError("No client connected")

        # Check client status
        status = self._get_status()

        if status == MSG_NEEDINIT:
            self._send_init(atoms)
            status = self._get_status()

        if status != MSG_READY:
            raise RuntimeError(f"Unexpected status: {status!r}")

        # Send positions
        self._send_posdata(atoms)

        # Wait for calculation
        status = self._get_status()
        if status != MSG_HAVEDATA:
            raise RuntimeError(f"Expected HAVEDATA, got: {status!r}")

        # Get forces
        return self._get_forces(len(atoms))

    def close(self) -> None:
        """Close the server and client connections."""
        if self._client:
            try:
                self._send_message(MSG_EXIT)
            except Exception:
                pass
            self._client.close()
            self._client = None

        if self._socket:
            self._socket.close()
            self._socket = None

    def _send_message(self, msg: bytes) -> None:
        """Send a protocol message."""
        if len(msg) != HEADER_SIZE:
            msg = msg.ljust(HEADER_SIZE)[:HEADER_SIZE]
        self._client.sendall(msg)

    def _recv_message(self) -> bytes:
        """Receive a protocol message."""
        data = b""
        while len(data) < HEADER_SIZE:
            chunk = self._client.recv(HEADER_SIZE - len(data))
            if not chunk:
                raise ConnectionError("Connection closed")
            data += chunk
        return data

    def _get_status(self) -> bytes:
        """Request and receive client status."""
        self._send_message(MSG_STATUS)
        return self._recv_message()

    def _send_init(self, atoms: Atoms) -> None:
        """Send initialization data."""
        self._send_message(MSG_INIT)

        # Send bead index (0 for single calculation)
        self._client.sendall(struct.pack("i", 0))

        # Send cell (not used, but protocol requires it)
        init_data = b"init"
        self._client.sendall(struct.pack("i", len(init_data)))
        self._client.sendall(init_data)

        self._initialized = True

    def _send_posdata(self, atoms: Atoms) -> None:
        """Send cell vectors and atomic positions."""
        self._send_message(MSG_POSDATA)

        # Cell vectors (in Bohr, row-major)
        cell = atoms.get_cell() * ANGSTROM_TO_BOHR
        cell_flat = cell.T.flatten()  # Column-major for i-PI
        self._client.sendall(cell_flat.astype(">f8").tobytes())

        # Inverse cell
        inv_cell = np.linalg.inv(cell.T).flatten()
        self._client.sendall(inv_cell.astype(">f8").tobytes())

        # Number of atoms
        n_atoms = len(atoms)
        self._client.sendall(struct.pack(">i", n_atoms))

        # Positions (in Bohr)
        positions = atoms.get_positions().flatten() * ANGSTROM_TO_BOHR
        self._client.sendall(positions.astype(">f8").tobytes())

    def _get_forces(self, n_atoms: int) -> dict:
        """Request and receive forces."""
        self._send_message(MSG_GETFORCE)

        # Wait for FORCEREADY
        status = self._recv_message()
        if status != MSG_FORCEREADY:
            raise RuntimeError(f"Expected FORCEREADY, got: {status!r}")

        # Energy (Hartree -> eV)
        energy_data = self._client.recv(8)
        energy = struct.unpack(">d", energy_data)[0] * HARTREE_TO_EV

        # Number of atoms
        natoms_data = self._client.recv(4)
        n_recv = struct.unpack(">i", natoms_data)[0]
        if n_recv != n_atoms:
            raise RuntimeError(f"Atom count mismatch: {n_recv} vs {n_atoms}")

        # Forces (Hartree/Bohr -> eV/Angstrom)
        forces_data = self._client.recv(n_atoms * 3 * 8)
        forces = np.frombuffer(forces_data, dtype=">f8").reshape(n_atoms, 3)
        forces = forces * HARTREE_TO_EV / BOHR_TO_ANGSTROM

        # Virial/stress (9 components)
        virial_data = self._client.recv(72)
        virial = np.frombuffer(virial_data, dtype=">f8").reshape(3, 3)

        # Extra bytes (i-PI protocol)
        extra_len_data = self._client.recv(4)
        extra_len = struct.unpack(">i", extra_len_data)[0]
        if extra_len > 0:
            self._client.recv(extra_len)

        return {
            "energy": energy,
            "forces": forces,
            "virial": virial,
        }

    def __enter__(self) -> SocketServer:
        self.start()
        return self

    def __exit__(self, *args) -> None:
        self.close()


# =============================================================================
# Socket Client (VASP Side)
# =============================================================================


class SocketClient:
    """Socket client for VASP to receive positions and send forces.

    This is typically run by a wrapper around VASP. It receives
    positions from a driver, passes them to VASP, and sends
    back energies and forces.

    Args:
        config: Socket configuration.
        runner: InteractiveRunner instance for VASP communication.

    Example:
        >>> from vasp.runners import InteractiveRunner
        >>>
        >>> runner = InteractiveRunner()
        >>> client = SocketClient(config, runner)
        >>> client.connect()
        >>> client.run()  # Main loop
    """

    def __init__(
        self,
        config: SocketConfig,
        runner: InteractiveRunner | None = None,
    ):
        self.config = config
        self.runner = runner
        self._socket: socket.socket | None = None
        self._state = SocketState.DISCONNECTED
        self._atoms: Atoms | None = None
        self._results: dict | None = None

    def connect(self) -> None:
        """Connect to the server."""
        if self.config.unix_socket:
            self._socket = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            self._socket.connect(self.config.unix_socket)
        else:
            self._socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self._socket.connect((self.config.host, self.config.port))

        self._socket.settimeout(self.config.timeout)
        self._state = SocketState.CONNECTED

    def run(self, atoms_template: Atoms, directory: str) -> None:
        """Main client loop.

        Args:
            atoms_template: Template atoms object (for element types).
            directory: VASP calculation directory.
        """

        self._atoms = atoms_template.copy()
        first_step = True

        while True:
            msg = self._recv_message()

            if msg == MSG_STATUS:
                if self._results is not None:
                    self._send_message(MSG_HAVEDATA)
                elif first_step:
                    self._send_message(MSG_NEEDINIT)
                else:
                    self._send_message(MSG_READY)

            elif msg == MSG_INIT:
                # Receive init data (bead index + init string)
                struct.unpack("i", self._socket.recv(4))[0]
                init_len = struct.unpack("i", self._socket.recv(4))[0]
                if init_len > 0:
                    self._socket.recv(init_len)
                first_step = False

            elif msg == MSG_POSDATA:
                # Receive cell and positions
                cell, positions = self._recv_posdata(len(self._atoms))

                # Update atoms
                self._atoms.set_cell(cell)
                self._atoms.set_positions(positions)
                self._atoms.set_pbc(True)

                # Run VASP
                if first_step:
                    results = self.runner.start(directory, self._atoms)
                    first_step = False
                else:
                    results = self.runner.step(self._atoms)

                # Store results for GETFORCE
                self._results = {
                    "energy": results.energy,
                    "forces": results.forces,
                    "stress": results.stress,
                }

            elif msg == MSG_GETFORCE:
                if self._results is None:
                    raise RuntimeError("GETFORCE without calculation")
                self._send_forces()
                self._results = None

            elif msg == MSG_EXIT:
                break

            else:
                raise RuntimeError(f"Unknown message: {msg!r}")

        self.close()

    def close(self) -> None:
        """Close the connection."""
        if self._socket:
            self._socket.close()
            self._socket = None
        self._state = SocketState.DISCONNECTED

        if self.runner:
            self.runner.close()

    def _send_message(self, msg: bytes) -> None:
        """Send a protocol message."""
        if len(msg) != HEADER_SIZE:
            msg = msg.ljust(HEADER_SIZE)[:HEADER_SIZE]
        self._socket.sendall(msg)

    def _recv_message(self) -> bytes:
        """Receive a protocol message."""
        data = b""
        while len(data) < HEADER_SIZE:
            chunk = self._socket.recv(HEADER_SIZE - len(data))
            if not chunk:
                raise ConnectionError("Connection closed")
            data += chunk
        return data

    def _recv_posdata(self, n_atoms: int) -> tuple[np.ndarray, np.ndarray]:
        """Receive cell and positions."""
        # Cell (9 doubles, Bohr -> Angstrom)
        cell_data = self._socket.recv(72)
        cell = np.frombuffer(cell_data, dtype=">f8").reshape(3, 3).T
        cell = cell * BOHR_TO_ANGSTROM

        # Inverse cell (ignored)
        self._socket.recv(72)

        # Number of atoms
        n_data = self._socket.recv(4)
        n_recv = struct.unpack(">i", n_data)[0]
        if n_recv != n_atoms:
            raise RuntimeError(f"Atom count mismatch: {n_recv} vs {n_atoms}")

        # Positions (Bohr -> Angstrom)
        pos_data = self._socket.recv(n_atoms * 3 * 8)
        positions = np.frombuffer(pos_data, dtype=">f8").reshape(n_atoms, 3)
        positions = positions * BOHR_TO_ANGSTROM

        return cell, positions

    def _send_forces(self) -> None:
        """Send energy, forces, and virial."""
        self._send_message(MSG_FORCEREADY)

        assert self._atoms is not None
        assert self._results is not None

        n_atoms = len(self._atoms)

        # Energy (eV -> Hartree)
        energy = self._results["energy"] * EV_TO_HARTREE
        self._socket.sendall(struct.pack(">d", energy))

        # Number of atoms
        self._socket.sendall(struct.pack(">i", n_atoms))

        # Forces (eV/Angstrom -> Hartree/Bohr)
        forces = self._results["forces"] * EV_TO_HARTREE * BOHR_TO_ANGSTROM
        self._socket.sendall(forces.astype(">f8").tobytes())

        # Virial (convert from stress if available)
        virial = np.zeros((3, 3))
        if self._results.get("stress") is not None:
            # VASP stress is in eV/Å³ (Voigt: xx, yy, zz, yz, xz, xy)
            # Virial = -stress * volume, in Hartree
            stress = self._results["stress"]
            volume = self._atoms.get_volume()  # Å³

            # Convert Voigt to 3x3 tensor
            # Voigt order: xx, yy, zz, yz, xz, xy
            virial[0, 0] = -stress[0] * volume  # xx
            virial[1, 1] = -stress[1] * volume  # yy
            virial[2, 2] = -stress[2] * volume  # zz
            virial[1, 2] = virial[2, 1] = -stress[3] * volume  # yz
            virial[0, 2] = virial[2, 0] = -stress[4] * volume  # xz
            virial[0, 1] = virial[1, 0] = -stress[5] * volume  # xy

            # Convert eV -> Hartree
            virial *= EV_TO_HARTREE
        self._socket.sendall(virial.astype(">f8").tobytes())

        # Extra bytes (empty)
        self._socket.sendall(struct.pack(">i", 0))

    def __enter__(self) -> SocketClient:
        self.connect()
        return self

    def __exit__(self, *args) -> None:
        self.close()
