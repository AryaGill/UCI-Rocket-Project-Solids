import time
from PyQt6.QtCore import QThread, pyqtSignal


def list_serial_ports():
    """Return a list of (device, description) tuples for all available serial ports."""
    try:
        from serial.tools import list_ports
    except ImportError:
        print("ERROR: pyserial not installed. Run: pip install pyserial")
        return []

    ports = []
    for p in list_ports.comports():
        ports.append((p.device, p.description))
    return ports


def select_serial_port():
    """
    Print all available serial ports and prompt user to select one.
    Returns: Selected port device string or None if skipped/no ports.
    """
    ports = list_serial_ports()

    if not ports:
        print("No serial ports found.")
        print("Make sure your device is connected and drivers are installed.")
        return None

    print("\n" + "="*60)
    print("Available Serial Ports:")
    print("="*60)

    for i, (device, description) in enumerate(ports):
        print(f"  [{i}] {device}")
        print(f"      {description}")

    print("="*60)

    while True:
        try:
            selection = input("\nSelect port number (or press Enter to skip): ").strip()

            if selection == "":
                print("Skipping serial connection.\n")
                return None

            idx = int(selection)
            if 0 <= idx < len(ports):
                selected = ports[idx][0]
                print(f"Selected: {selected}\n")
                return selected
            else:
                print(f"Invalid selection. Please enter a number between 0 and {len(ports)-1}.")

        except ValueError:
            print("Invalid input. Please enter a number or press Enter to skip.")
        except KeyboardInterrupt:
            print("\nCancelled by user.")
            return None


class SerialStreamer(QThread):
    """
    Thread that connects to a serial port, reads CSV telemetry data line-by-line,
    and emits data dictionaries to update frontend graphs.

    Expected CSV format (no header required):
    Time, Temp, Pressure, Alt, Gyro_X, Gyro_Y, Gyro_Z, Accel_X1, Accel_Y1, Accel_Z1, Accel_X2, Accel_Y2, Accel_Z2, flight_state

    Example line:
    1.523, 25.3, 101325.0, 123.5, 0.01, -0.02, 0.03, 0.98, 0.02, 9.81, 1.2, -0.5, 15.3, 2
    """

    # Fixed column mapping - position-based, not header-based
    COLUMNS = [
        "Time",
        "Temp",
        "Pressure",
        "Alt",
        "Filtered_Alt",
        "Gyro_X",
        "Gyro_Y",
        "Gyro_Z",
        "Accel_X1",
        "Accel_Y1",
        "Accel_Z1",
        "Accel_world_x",
        "Accel_world_y",
        "Accel_world_z",
        "mag_r",
        "mag_p",
        "mag_y",
        "roll",
        "pitch",
        "yaw",
        "velocity_x",
        "velocity_y",
        "velocity_z",
        "Quaternion_W",
        "Quaternion_X",
        "Quaternion_Y",
        "Quaternion_Z",
        "pred_apo",
        "AB_Deployment",
        "Cam1V",
        "Cam2V",
        "main_p_ematch_voltage",
        "main_s_ematch_voltage",
        "drogue_p_ematch_voltage",
        "drogue_s_ematch_voltage",
        "flight_state",
        "command_echo"
    ]
    
    new_data = pyqtSignal(dict)
    finished = pyqtSignal()
    status = pyqtSignal(str)
    fc_message = pyqtSignal(str)

    def __init__(self, port, baud=57600, timeout=0.5, parent=None):
        """
        Initialize serial streamer.

        Args:
            port: Serial port device (e.g., '/dev/ttyUSB0' or 'COM3')
            baud: Baud rate (default: 115200)
            timeout: Read timeout in seconds (default: 0.5)
            parent: Parent QObject
        """
        super().__init__(parent)
        self.port = port
        self.baud = baud
        self.timeout = timeout
        self.is_running = True
        self.paused = False
        self._ser = None

    def run(self):
        """Main thread loop - opens serial connection and reads data continuously."""
        print(f"Starting SerialStreamer thread for port: {self.port}")
        try:
            import serial

            self.status.emit(f"Connecting to {self.port} @ {self.baud} baud...")

            self._ser = serial.Serial(
                port=self.port,
                baudrate=self.baud,
                timeout=self.timeout,
                bytesize=serial.EIGHTBITS,
                parity=serial.PARITY_NONE,
                stopbits=serial.STOPBITS_ONE
            )

            time.sleep(0.5)
            self._ser.reset_input_buffer()

            self.status.emit(f"✓ Connected to {self.port}")
            self.status.emit(f"Expecting {len(self.COLUMNS)} columns: {', '.join(self.COLUMNS)}")

            while self.is_running:
                while self.paused and self.is_running:
                    time.sleep(0.05)

                if not self.is_running:
                    break

                try:
                    line = self._ser.readline()
                except Exception as e:
                    self.status.emit(f"Read error: {e}")
                    continue

                if not line:
                    continue

                # print(f"RAW BYTES: {line!r}") 


                try:
                    text = line.decode('utf-8', errors='ignore').strip()
                except Exception:
                    continue

                if not text:
                    continue

                data = self._parse_line(text)
                if data:
                    # print(f"Parsed data: {data}")
                    self.new_data.emit(data)

        except Exception as e:
            self.status.emit(f"Serial connection error: {e}")
        finally:
            if self._ser:
                try:
                    self._ser.close()
                    self.status.emit(f"Disconnected from {self.port}")
                except Exception:
                    pass
            self.finished.emit()

    def _parse_line(self, text: str) -> dict | None:
        """
        Parse a CSV line into a dictionary using fixed column positions.

        Expected format:
        Time, Temp, Pressure, Alt, Gyro_X, Gyro_Y, Gyro_Z, Accel_X1, Accel_Y1, Accel_Z1, Accel_X2, Accel_Y2, Accel_Z2,  flight_state

        Returns:
            Dictionary mapping column names to values, or None if invalid
        """
        parts = [p.strip() for p in text.split(",")]

        if len(parts) != len(self.COLUMNS):
            # Treat it as a plain-text FC message (ack, confirmation, etc.)
            self.fc_message.emit(text)
            return None

        data = {}
        for col_name, value_str in zip(self.COLUMNS, parts):
            if col_name != "command_echo":
                data[col_name] = self._to_number(value_str)
            else:
                data[col_name] = value_str

        return data

    @staticmethod
    def _to_number(s: str):
        """Convert string to int/float if possible, otherwise return string."""
        try:
            if "." in s or "e" in s.lower():
                return float(s)
            return int(s)
        except Exception:
            return s

    #The whole write command was rewritten using non-blocking
    def write_command(self, cmd: str, burst: int = 20,
                  tx_sleep: float = 0.05, rx_window: float = 0.20):
        """
        Send a command over a half-duplex LoRa radio using a strict
        alternating TX / RX pattern:

            For each of `burst` repetitions:
                1. Pause the RX read-loop
                2. Transmit one packet
                3. Resume the RX read-loop for `rx_window` seconds
                so any incoming telemetry is not missed
                4. Repeat

        Args:
            cmd:       Command string (newline appended automatically).
            burst:     Number of times to repeat the command (default 20).
            tx_sleep:  Seconds to wait after writing before opening the RX
                    window — gives the radio time to finish transmitting.
            rx_window: Seconds to leave the read-loop running between
                    transmissions so telemetry packets are not dropped.
        """
        if not (self._ser and self._ser.is_open):
            self.status.emit("Write error: port not open")
            return

        encoded = (cmd.strip() + "\n").encode("utf-8")

        def _burst():
            was_paused = self.paused          # preserve caller's pause state

            try:
                for i in range(burst):
                    if not self.is_running or not self._ser.is_open:
                        break

                    # ── 1. PAUSE RX ──────────────────────────────────────────
                    # Signal the read loop to stop, then wait long enough for
                    # any blocking readline() that is already in progress to
                    # return (timeout is 0.5 s by default).
                    self.paused = True
                    time.sleep(self.timeout + 0.02)   # outlast the current readline()

                    # ── 2. TRANSMIT ──────────────────────────────────────────
                    try:
                        self._ser.write(encoded)
                        self.status.emit(f"TX [{i+1}/{burst}]: {cmd.strip()}")
                    except Exception as e:
                        self.status.emit(f"Write error on burst {i+1}: {e}")
                        break

                    time.sleep(tx_sleep)   # let the radio finish transmitting

                    # ── 3. OPEN RX WINDOW ────────────────────────────────────
                    # Resume the read-loop so the main run() can process any
                    # telemetry the rocket sends back between our transmissions.
                    self.paused = False
                    time.sleep(rx_window)  # receive window before next TX

                self.status.emit(f"Burst complete – sent {burst}× '{cmd.strip()}'")

            finally:
                # Always restore the pause state that was active before the
                # command was issued (e.g. if the user had manually paused).
                self.paused = was_paused

        import threading
        threading.Thread(target=_burst, daemon=True).start()

    def stop(self):
        """Stop the streaming thread."""
        self.is_running = False

    def pause(self):
        """Pause data streaming (connection remains open)."""
        self.paused = True
        self.status.emit("Paused")

    def resume(self):
        """Resume data streaming."""
        self.paused = False
        self.status.emit("Resumed")


if __name__ == "__main__":
    # Test the port selection
    port = select_serial_port()
    if port:
        print(f"Would connect to: {port}")
        print(f"\nExpecting CSV format with {len(SerialStreamer.COLUMNS)} columns:")
        for i, col in enumerate(SerialStreamer.COLUMNS, 1):
            print(f"  Column {i}: {col}")