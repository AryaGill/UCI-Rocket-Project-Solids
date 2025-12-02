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
