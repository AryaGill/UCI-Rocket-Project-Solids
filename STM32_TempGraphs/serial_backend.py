"""
Serial communication backend
Handles serial port discovery, connection, and data reading
"""

import serial
import serial.tools.list_ports
import csv
from io import StringIO
from sensor_config import SENSORS


class SerialBackend:
    def __init__(self, sensor_type='LPS22HHTR', baudrate=115200):
        self.sensor_type = sensor_type
        self.baudrate = baudrate
        self.serial_conn = None
        self.expected_columns = SENSORS[sensor_type]['columns']

    @staticmethod
    def list_serial_ports():
        """List all available serial ports"""
        ports = serial.tools.list_ports.comports()
        available_ports = []

        print("\nAvailable Serial Ports:")
        print("-" * 50)

        for idx, port in enumerate(ports, 1):
            print(f"{idx}. {port.device}")
            print(f"   Description: {port.description}")
            print(f"   Hardware ID: {port.hwid}")
            print()
            available_ports.append(port.device)

        if not available_ports:
            print("No serial ports found!")

        return available_ports

    def connect(self, port, baudrate=None):
        """Connect to specified serial port"""
        if baudrate is None:
            baudrate = self.baudrate

        try:
            self.serial_conn = serial.Serial(
                port=port,
                baudrate=baudrate,
                timeout=1,
                bytesize=serial.EIGHTBITS,
                parity=serial.PARITY_NONE,
                stopbits=serial.STOPBITS_ONE
            )
            print(f"\n✓ Connected to {port} at {baudrate} baud")
            return True
        except serial.SerialException as e:
            print(f"\n✗ Failed to connect to {port}: {e}")
            return False

    def validate_data(self, data_row):
        """Validate that received data matches expected format"""
        if len(data_row) != len(self.expected_columns):
            return False, f"Expected {len(self.expected_columns)} columns, got {len(data_row)}"

        # Check if all values are numeric
        try:
            [float(val) for val in data_row]
            return True, None
        except ValueError as e:
            return False, f"Non-numeric data received: {e}"

    def read_line(self):
        """Read and parse one line of CSV data from serial"""
        if not self.serial_conn or not self.serial_conn.is_open:
            return None, "Serial connection not open"

        try:
            line = self.serial_conn.readline().decode('utf-8').strip()

            if not line:
                return None, "Empty line"

            # Parse CSV
            reader = csv.reader(StringIO(line))
            data_row = next(reader)

            # Validate
            is_valid, error_msg = self.validate_data(data_row)

            if is_valid:
                # Convert to float and create dict
                data_dict = {
                    col: float(val) 
                    for col, val in zip(self.expected_columns, data_row)
                }
                return data_dict, None
            else:
                return None, error_msg

        except Exception as e:
            return None, f"Error reading serial: {e}"

    def disconnect(self):
        """Close serial connection"""
        if self.serial_conn and self.serial_conn.is_open:
            self.serial_conn.close()
            print("\n✓ Serial connection closed")

    def __del__(self):
        self.disconnect()


def select_serial_port():
    """Interactive prompt to select a serial port"""
    backend = SerialBackend()
    ports = backend.list_serial_ports()

    if not ports:
        return None

    while True:
        try:
            choice = input(f"Select port (1-{len(ports)}): ").strip()
            idx = int(choice) - 1

            if 0 <= idx < len(ports):
                return ports[idx]
            else:
                print(f"Please enter a number between 1 and {len(ports)}")
        except ValueError:
            print("Invalid input. Please enter a number.")
        except KeyboardInterrupt:
            print("\nCancelled")
            return None


if __name__ == "__main__":
    # Test the backend
    port = select_serial_port()

    if port:
        backend = SerialBackend(sensor_type='LPS22HHTR')

        if backend.connect(port):
            print("\nReading data (Press Ctrl+C to stop)...")
            print(f"Expected format: {', '.join(backend.expected_columns)}")
            print("-" * 50)

            try:
                while True:
                    data, error = backend.read_line()

                    if data:
                        print(data)
                    elif error and "Empty line" not in error:
                        print(f"Error: {error}")

            except KeyboardInterrupt:
                print("\n\nStopped by user")
            finally:
                backend.disconnect()