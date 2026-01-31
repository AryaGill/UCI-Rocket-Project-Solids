"""
Serial Port Monitor
Lists available serial ports, connects to selected port, and prints incoming data.
"""

import serial
import serial.tools.list_ports
import sys

def list_serial_ports():
    """List all available serial ports."""
    ports = serial.tools.list_ports.comports()
    return ports

def main():
    print("=== Serial Port Monitor ===\n")

    # List available ports
    ports = list_serial_ports()

    if not ports:
        print("No serial ports found!")
        sys.exit(1)

    print("Available serial ports:")
    for i, port in enumerate(ports, 1):
        print(f"{i}. {port.device} - {port.description}")

    # Get user selection
    while True:
        try:
            choice = input(f"\nSelect port (1-{len(ports)}): ")
            port_index = int(choice) - 1

            if 0 <= port_index < len(ports):
                selected_port = ports[port_index].device
                break
            else:
                print(f"Please enter a number between 1 and {len(ports)}")
        except ValueError:
            print("Please enter a valid number")
        except KeyboardInterrupt:
            print("\nExiting...")
            sys.exit(0)

    # Get baud rate
    try:
        baud_rate = input("\nEnter baud rate (default 9600): ").strip()
        baud_rate = int(baud_rate) if baud_rate else 9600
    except ValueError:
        print("Invalid baud rate, using 9600")
        baud_rate = 9600

    # Connect to selected port
    print(f"\nConnecting to {selected_port} at {baud_rate} baud...")

    try:
        ser = serial.Serial(
            port=selected_port,
            baudrate=baud_rate,
            timeout=1,
            bytesize=serial.EIGHTBITS,
            parity=serial.PARITY_NONE,
            stopbits=serial.STOPBITS_ONE
        )

        print(f"Connected to {selected_port}")
        print("Reading data (Press Ctrl+C to stop)...\n")
        print("-" * 50)

        # Read and print data continuously
        while True:
            if ser.in_waiting > 0:
                try:
                    data = ser.readline().decode('utf-8', errors='replace').rstrip()
                    if data:
                        print(data)
                except UnicodeDecodeError:
                    # If decode fails, print raw bytes
                    raw_data = ser.readline()
                    print(f"[RAW] {raw_data}")

    except serial.SerialException as e:
        print(f"Error: Could not open port {selected_port}")
        print(f"Details: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        print("\n" + "-" * 50)
        print("\nStopped monitoring. Closing connection...")
        ser.close()
        print("Connection closed.")
        sys.exit(0)
    except Exception as e:
        print(f"\nUnexpected error: {e}")
        if 'ser' in locals() and ser.is_open:
            ser.close()
        sys.exit(1)

if __name__ == "__main__":
    main()