import pandas as pd
import time
from PyQt6.QtCore import QThread, pyqtSignal


class DataStreamer(QThread):
    """
    Backend thread that reads CSV file row-by-row and emits data with delays
    to simulate real-time rocket telemetry.
    """
    new_data = pyqtSignal(dict)  # Signal emits dictionary of column:value pairs
    finished = pyqtSignal()
    
    def __init__(self, csv_file, delay=0.1):
        """
        Initialize the data streamer.
        
        Args:
            csv_file (str): Path to CSV file with telemetry data
            delay (float): Delay between rows in seconds (default 0.1s = 100ms)
        """
        super().__init__()
        self.csv_file = csv_file
        self.delay = delay
        self.is_running = True
        self.paused = False
        
    def run(self):
        """Read CSV file and emit data row by row with delays."""
        try:
            # Read CSV file
            df = pd.read_csv(self.csv_file)
            
            # Iterate through each row
            for index, row in df.iterrows():
                if not self.is_running:
                    break
                    
                # Wait if paused
                while self.paused and self.is_running:
                    time.sleep(0.1)
                    
                if not self.is_running:
                    break
                
                # Convert row to dictionary and emit
                data_dict = row.to_dict()
                self.new_data.emit(data_dict)
                
                # Delay to simulate real-time data
                time.sleep(self.delay)
            
            self.finished.emit()
            
        except Exception as e:
            print(f"Error reading CSV: {e}")
            self.finished.emit()
    
    def stop(self):
        """Stop the data streaming."""
        self.is_running = False
        
    def pause(self):
        """Pause the data streaming."""
        self.paused = True
        
    def resume(self):
        """Resume the data streaming."""
        self.paused = False

def list_serial_ports():
    """Return a list of (device, description) tuples."""
    try:
        from serial.tools import list_ports
    except Exception:
        return []
    ports = []
    for p in list_ports.comports():
        ports.append((p.device, p.description))
    return ports


class SerialStreamer(QThread):
    """
    Reads newline-delimited telemetry from a serial port and emits dicts.
    Also provides a thread-safe-ish write() to send commands back.
    """
    new_data = pyqtSignal(dict)
    finished = pyqtSignal()
    status = pyqtSignal(str)

    def __init__(self, port: str, baud: int = 115200, timeout: float = 0.2, parent=None):
        super().__init__(parent)
        self.port = port
        self.baud = baud
        self.timeout = timeout
        self.is_running = True
        self.paused = False
        self._ser = None

        # If telemetry includes a header line, we’ll store it here.
        self._columns = None

    def run(self):
        try:
            import serial
            self._ser = serial.Serial(self.port, self.baud, timeout=self.timeout)
            self.status.emit(f"Connected to {self.port} @ {self.baud}")

            while self.is_running:
                # pause behavior consistent with your CSV streamer
                while self.paused and self.is_running:
                    time.sleep(0.05)
                if not self.is_running:
                    break

                line = self._ser.readline()
                if not line:
                    continue

                try:
                    text = line.decode(errors="ignore").strip()
                except Exception:
                    continue

                if not text:
                    continue

                data = self._parse_line(text)
                if data:
                    self.new_data.emit(data)

        except Exception as e:
            self.status.emit(f"Serial error: {e}")
        finally:
            try:
                if self._ser:
                    self._ser.close()
            except Exception:
                pass
            self.finished.emit()

    def _parse_line(self, text: str) -> dict | None:
        """
        Supports either:
        1) CSV with header first: "Time,Alt,Temp,MagX,..."
           then data:          "0.1,123,24.0, ..."

        2) Key=Value pairs: "Time=0.1 Alt=123 Temp=24.0 ..."
        """
        # Key=Value format
        if "=" in text and ("," not in text):
            out = {}
            parts = text.split()
            for p in parts:
                if "=" not in p:
                    continue
                k, v = p.split("=", 1)
                out[k.strip()] = self._to_number(v.strip())
            return out if out else None

        # CSV format
        parts = [p.strip() for p in text.split(",")]
        if not parts:
            return None

        # detect header (non-numeric tokens)
        if self._columns is None:
            looks_like_header = any(self._to_number(x) is None for x in parts)
            if looks_like_header:
                self._columns = parts
                self.status.emit("Telemetry header received")
                return None
            # if no header, you can hard-map columns here if needed
            return None  # force header for now

        if len(parts) != len(self._columns):
            return None

        return {k: self._to_number(v) for k, v in zip(self._columns, parts)}

    @staticmethod
    def _to_number(s: str):
        try:
            if "." in s or "e" in s.lower():
                return float(s)
            return int(s)
        except Exception:
            return s

    def write_command(self, cmd: str):
        """Send a line back out over the same serial link (USB or RF modem)."""
        try:
            if self._ser and self._ser.is_open:
                self._ser.write((cmd.strip() + "\n").encode())
        except Exception:
            pass

    def stop(self):
        self.is_running = False

    def pause(self):
        self.paused = True

    def resume(self):
        self.paused = False

