"""
Serial Test Data Generator - Simulates STM32 serial output
Creates a virtual serial connection for testing GUI without hardware
"""

import serial
import time
import math
import random
import sys
import argparse


def generate_imu_data():
    """Generate realistic ICM45686 IMU test data"""
    t = time.time()

    # Simulate some motion patterns
    accel_x = math.sin(t) * 2.0 + random.gauss(0, 0.1)
    accel_y = math.cos(t) * 1.5 + random.gauss(0, 0.1)
    accel_z = -9.81 + random.gauss(0, 0.2)  # Gravity with noise

    gyro_x = math.sin(t * 2) * 0.5 + random.gauss(0, 0.05)
    gyro_y = math.cos(t * 2) * 0.3 + random.gauss(0, 0.05)
    gyro_z = random.gauss(0, 0.02)

    # Format as CSV
    csv_line = f"{accel_x:.6f},{accel_y:.6f},{accel_z:.6f},{gyro_x:.6f},{gyro_y:.6f},{gyro_z:.6f}\n"
    return csv_line


def run_serial_simulator(port, baudrate=115200, rate=10):
    """
    Send test data to a serial port

    Args:
        port: Serial port name (e.g., 'COM3', '/dev/ttyUSB0', or '/dev/pts/X')
        baudrate: Baud rate (default 115200)
        rate: Data rate in Hz (default 10)
    """
    try:
        ser = serial.Serial(port, baudrate, timeout=1)
        print(f"✓ Connected to {port} at {baudrate} baud")
        print(f"✓ Sending data at {rate} Hz")
        print(f"✓ Press Ctrl+C to stop\n")
        print("-" * 60)

        counter = 0
        while True:
            data = generate_imu_data()
            ser.write(data.encode('utf-8'))

            # Print to console every 10 samples
            if counter % 10 == 0:
                print(f"Sent: {data.strip()}")

            counter += 1
            time.sleep(1.0 / rate)

    except serial.SerialException as e:
        print(f"✗ Serial error: {e}")
        print("\nMake sure the virtual serial port pair is created.")
        print("See instructions below.")
        sys.exit(1)
    except KeyboardInterrupt:
        print("\n\n✓ Stopped by user")
        if ser.is_open:
            ser.close()
            print("✓ Serial port closed")
    except Exception as e:
        print(f"✗ Error: {e}")
        sys.exit(1)


def print_usage_instructions():
    """Print instructions for setting up virtual serial ports"""
    print("""
╔════════════════════════════════════════════════════════════════╗
║          Virtual Serial Port Setup Instructions               ║
╔════════════════════════════════════════════════════════════════╗

This script sends test data to a serial port to simulate your STM32.
You need to create a VIRTUAL SERIAL PORT PAIR to test without hardware.

═══════════════════════════════════════════════════════════════════
OPTION 1: macOS/Linux - Using socat (Recommended)
═══════════════════════════════════════════════════════════════════

1. Install socat:
   macOS:  brew install socat
   Linux:  sudo apt-get install socat

2. Create virtual serial port pair:
   socat -d -d pty,raw,echo=0 pty,raw,echo=0

   This will output something like:
   2024/01/17 13:20:00 socat[12345] N PTY is /dev/ttys001
   2024/01/17 13:20:00 socat[12345] N PTY is /dev/ttys002

3. In one terminal, run this script with the FIRST port:
   python serial_test_generator.py /dev/ttys001

4. In another terminal/window, run the GUI and connect to SECOND port:
   python main_gui.py
   (Select /dev/ttys002 in the GUI)

═══════════════════════════════════════════════════════════════════
OPTION 2: Windows - Using com0com
═══════════════════════════════════════════════════════════════════

1. Download and install com0com:
   https://sourceforge.net/projects/com0com/

2. After installation, it creates COM port pairs (e.g., COM3 <-> COM4)

3. Run this script with one port:
   python serial_test_generator.py COM3

4. Run the GUI and connect to the paired port (COM4)

═══════════════════════════════════════════════════════════════════
OPTION 3: Linux - Using built-in pts
═══════════════════════════════════════════════════════════════════

Use socat as shown in OPTION 1 (works best on Linux)

═══════════════════════════════════════════════════════════════════
OPTION 4: Test with stdout (no GUI test)
═══════════════════════════════════════════════════════════════════

Just see the data being generated:
   python serial_test_generator.py --stdout

═══════════════════════════════════════════════════════════════════
""")


def main():
    parser = argparse.ArgumentParser(
        description='ICM45686 Serial Test Data Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        'port',
        nargs='?',
        help='Serial port name (e.g., COM3, /dev/ttyUSB0, /dev/ttys001)'
    )

    parser.add_argument(
        '--baudrate', '-b',
        type=int,
        default=115200,
        help='Baud rate (default: 115200)'
    )

    parser.add_argument(
        '--rate', '-r',
        type=int,
        default=10,
        help='Data rate in Hz (default: 10)'
    )

    parser.add_argument(
        '--stdout',
        action='store_true',
        help='Print to stdout instead of serial port'
    )

    parser.add_argument(
        '--help-setup',
        action='store_true',
        help='Show virtual serial port setup instructions'
    )

    args = parser.parse_args()

    if args.help_setup:
        print_usage_instructions()
        sys.exit(0)

    if args.stdout:
        print("Generating test data to stdout (Ctrl+C to stop)...")
        print(f"Rate: {args.rate} Hz")
        print("-" * 60)
        try:
            counter = 0
            while True:
                data = generate_imu_data()
                print(data.strip())
                counter += 1
                time.sleep(1.0 / args.rate)
        except KeyboardInterrupt:
            print(f"\n\nGenerated {counter} samples")
        sys.exit(0)

    if not args.port:
        print("✗ Error: No serial port specified\n")
        print("Usage: python serial_test_generator.py <port>")
        print("       python serial_test_generator.py --help-setup")
        print("       python serial_test_generator.py --stdout")
        print("\nExamples:")
        print("  python serial_test_generator.py /dev/ttys001")
        print("  python serial_test_generator.py COM3")
        print("  python serial_test_generator.py --stdout")
        sys.exit(1)

    run_serial_simulator(args.port, args.baudrate, args.rate)


if __name__ == "__main__":
    main()
