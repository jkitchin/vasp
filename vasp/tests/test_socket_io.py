"""Tests for socket I/O module."""

import socket
import struct
import threading

import numpy as np
import pytest

from vasp.runners.socket_io import (
    ANGSTROM_TO_BOHR,
    BOHR_TO_ANGSTROM,
    EV_TO_HARTREE,
    HARTREE_TO_EV,
    HEADER_SIZE,
    MSG_EXIT,
    MSG_FORCEREADY,
    MSG_GETFORCE,
    MSG_HAVEDATA,
    MSG_NEEDINIT,
    MSG_POSDATA,
    MSG_READY,
    MSG_STATUS,
    SocketClient,
    SocketConfig,
    SocketServer,
    SocketState,
)

# =============================================================================
# Unit Conversion Tests
# =============================================================================

class TestUnitConversions:
    """Test unit conversion constants."""

    def test_bohr_angstrom_roundtrip(self):
        """Test Bohr <-> Angstrom conversion roundtrip."""
        value = 5.43
        converted = value * ANGSTROM_TO_BOHR * BOHR_TO_ANGSTROM
        assert converted == pytest.approx(value)

    def test_hartree_ev_roundtrip(self):
        """Test Hartree <-> eV conversion roundtrip."""
        value = -10.5
        converted = value * EV_TO_HARTREE * HARTREE_TO_EV
        assert converted == pytest.approx(value)

    def test_bohr_value(self):
        """Test Bohr to Angstrom value."""
        # 1 Bohr ≈ 0.529 Angstrom
        assert BOHR_TO_ANGSTROM == pytest.approx(0.529, abs=0.001)

    def test_hartree_value(self):
        """Test Hartree to eV value."""
        # 1 Hartree ≈ 27.2 eV
        assert HARTREE_TO_EV == pytest.approx(27.2, abs=0.1)


# =============================================================================
# SocketConfig Tests
# =============================================================================

class TestSocketConfig:
    """Test SocketConfig dataclass."""

    def test_default_config(self):
        """Test default configuration values."""
        config = SocketConfig()

        assert config.host == "localhost"
        assert config.port == 31415
        assert config.unix_socket is None
        assert config.timeout == 3600.0

    def test_custom_config(self):
        """Test custom configuration."""
        config = SocketConfig(
            host="192.168.1.100",
            port=12345,
            timeout=60.0,
        )

        assert config.host == "192.168.1.100"
        assert config.port == 12345
        assert config.timeout == 60.0

    def test_unix_socket_config(self):
        """Test Unix socket configuration."""
        config = SocketConfig(unix_socket="/tmp/vasp.sock")

        assert config.unix_socket == "/tmp/vasp.sock"


# =============================================================================
# Protocol Message Tests
# =============================================================================

class TestProtocolMessages:
    """Test protocol message constants."""

    def test_message_lengths(self):
        """All messages should be 12 bytes."""
        messages = [
            MSG_STATUS, MSG_READY, MSG_NEEDINIT, MSG_HAVEDATA,
            MSG_POSDATA, MSG_GETFORCE, MSG_FORCEREADY, MSG_EXIT,
        ]

        for msg in messages:
            assert len(msg) == HEADER_SIZE

    def test_header_size(self):
        """Header size should be 12."""
        assert HEADER_SIZE == 12


# =============================================================================
# SocketState Tests
# =============================================================================

class TestSocketState:
    """Test SocketState enum."""

    def test_states(self):
        """Test all states exist."""
        assert SocketState.DISCONNECTED.value == "disconnected"
        assert SocketState.CONNECTED.value == "connected"
        assert SocketState.READY.value == "ready"
        assert SocketState.BUSY.value == "busy"


# =============================================================================
# SocketServer Tests
# =============================================================================

class TestSocketServer:
    """Test SocketServer class."""

    def test_init_default_config(self):
        """Test initialization with default config."""
        server = SocketServer()

        assert server.config.host == "localhost"
        assert server.config.port == 31415
        assert server._socket is None
        assert server._client is None

    def test_init_custom_config(self):
        """Test initialization with custom config."""
        config = SocketConfig(port=12345)
        server = SocketServer(config)

        assert server.config.port == 12345

    def test_start_creates_socket(self):
        """Test that start() creates a socket."""
        config = SocketConfig(port=0)  # Let OS choose port
        server = SocketServer(config)

        try:
            server.start()
            assert server._socket is not None
        finally:
            server.close()

    def test_close_cleans_up(self):
        """Test that close() cleans up resources."""
        config = SocketConfig(port=0)
        server = SocketServer(config)
        server.start()
        server.close()

        assert server._socket is None
        assert server._client is None

    def test_context_manager(self):
        """Test context manager protocol."""
        config = SocketConfig(port=0)

        with SocketServer(config) as server:
            assert server._socket is not None

        assert server._socket is None


# =============================================================================
# SocketClient Tests
# =============================================================================

class TestSocketClient:
    """Test SocketClient class."""

    def test_init(self):
        """Test client initialization."""
        config = SocketConfig()
        client = SocketClient(config)

        assert client.config == config
        assert client._socket is None
        assert client._state == SocketState.DISCONNECTED

    def test_initial_state(self):
        """Test initial state is disconnected."""
        config = SocketConfig()
        client = SocketClient(config)

        assert client._state == SocketState.DISCONNECTED


# =============================================================================
# Integration Tests
# =============================================================================

class TestServerClientIntegration:
    """Integration tests for server-client communication."""

    def test_connection(self):
        """Test basic client connection to server."""
        config = SocketConfig(port=0)  # Let OS choose port
        server = SocketServer(config)
        server.start()

        # Get actual port
        actual_port = server._socket.getsockname()[1]
        SocketConfig(port=actual_port)

        def connect_client():
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.connect(("localhost", actual_port))
            return sock

        # Start client in thread
        client_thread = threading.Thread(target=connect_client)
        client_thread.start()

        # Accept connection
        server.wait_for_client()

        assert server._client is not None

        client_thread.join(timeout=1)
        server.close()

    def test_message_exchange(self):
        """Test sending and receiving messages."""
        config = SocketConfig(port=0)
        server = SocketServer(config)
        server.start()
        actual_port = server._socket.getsockname()[1]

        responses = []

        def mock_client():
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.connect(("localhost", actual_port))

            # Receive STATUS
            msg = sock.recv(12)
            responses.append(msg)

            # Send READY
            sock.sendall(MSG_READY)

            # Receive EXIT
            msg = sock.recv(12)
            responses.append(msg)

            sock.close()

        client_thread = threading.Thread(target=mock_client)
        client_thread.start()

        server.wait_for_client()

        # Send STATUS
        server._send_message(MSG_STATUS)

        # Receive response
        response = server._recv_message()
        assert response == MSG_READY

        # Send EXIT
        server._send_message(MSG_EXIT)

        client_thread.join(timeout=2)
        server.close()

        assert MSG_STATUS in responses
        assert MSG_EXIT in responses


class TestBinaryProtocol:
    """Test binary data encoding/decoding."""

    def test_pack_int(self):
        """Test packing integers."""
        n = 42
        packed = struct.pack(">i", n)
        unpacked = struct.unpack(">i", packed)[0]

        assert unpacked == n

    def test_pack_double(self):
        """Test packing doubles."""
        value = -10.84274516
        packed = struct.pack(">d", value)
        unpacked = struct.unpack(">d", packed)[0]

        assert unpacked == pytest.approx(value)

    def test_pack_array(self):
        """Test packing numpy arrays."""
        arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        packed = arr.astype(">f8").tobytes()
        unpacked = np.frombuffer(packed, dtype=">f8")

        assert np.allclose(unpacked, arr)

    def test_cell_encoding(self):
        """Test cell vector encoding."""
        # 3x3 cell in Angstrom
        cell = np.array([
            [5.43, 0.0, 0.0],
            [0.0, 5.43, 0.0],
            [0.0, 0.0, 5.43],
        ])

        # Convert to Bohr and flatten (column-major for i-PI)
        cell_bohr = cell * ANGSTROM_TO_BOHR
        cell_flat = cell_bohr.T.flatten()

        # Pack and unpack
        packed = cell_flat.astype(">f8").tobytes()
        unpacked = np.frombuffer(packed, dtype=">f8").reshape(3, 3).T

        # Convert back to Angstrom
        recovered = unpacked * BOHR_TO_ANGSTROM

        assert np.allclose(recovered, cell)

    def test_forces_encoding(self):
        """Test forces encoding with unit conversion."""
        # Forces in eV/Angstrom
        forces = np.array([
            [0.1, -0.2, 0.3],
            [-0.1, 0.2, -0.3],
        ])

        # Convert to Hartree/Bohr
        forces_atomic = forces * EV_TO_HARTREE * BOHR_TO_ANGSTROM

        # Pack and unpack
        packed = forces_atomic.astype(">f8").tobytes()
        unpacked = np.frombuffer(packed, dtype=">f8").reshape(2, 3)

        # Convert back to eV/Angstrom
        recovered = unpacked * HARTREE_TO_EV / BOHR_TO_ANGSTROM

        assert np.allclose(recovered, forces)

    def test_stress_to_virial_conversion(self):
        """Test stress (Voigt) to virial (3x3) conversion."""
        # Stress in Voigt notation: xx, yy, zz, yz, xz, xy (eV/Å³)
        stress = np.array([0.1, 0.2, 0.3, 0.01, 0.02, 0.03])
        volume = 100.0  # Å³

        # Convert to virial: V_ij = -σ_ij * V
        virial = np.zeros((3, 3))
        virial[0, 0] = -stress[0] * volume  # xx
        virial[1, 1] = -stress[1] * volume  # yy
        virial[2, 2] = -stress[2] * volume  # zz
        virial[1, 2] = virial[2, 1] = -stress[3] * volume  # yz
        virial[0, 2] = virial[2, 0] = -stress[4] * volume  # xz
        virial[0, 1] = virial[1, 0] = -stress[5] * volume  # xy

        # Check symmetry
        assert np.allclose(virial, virial.T)

        # Check diagonal values
        assert virial[0, 0] == pytest.approx(-10.0)  # -0.1 * 100
        assert virial[1, 1] == pytest.approx(-20.0)  # -0.2 * 100
        assert virial[2, 2] == pytest.approx(-30.0)  # -0.3 * 100

        # Check off-diagonal
        assert virial[1, 2] == pytest.approx(-1.0)  # -0.01 * 100
        assert virial[0, 2] == pytest.approx(-2.0)  # -0.02 * 100
        assert virial[0, 1] == pytest.approx(-3.0)  # -0.03 * 100
