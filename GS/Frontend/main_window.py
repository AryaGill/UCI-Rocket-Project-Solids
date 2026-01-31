from PyQt6.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                             QGridLayout, QPushButton, QLabel)
from PyQt6.QtCore import Qt
from Backend.backend import SerialStreamer

class GroundStationWindow(QMainWindow):
    """
    Main window for the Ground Station application.
    Displays real-time telemetry graphs from the flight computer.
    """
    def __init__(self, port=None, parent=None):
        super().__init__(parent)
        self.selected_port = port
        self.streamer = None
        self.pyro_panel = None  # Will hold PyroPanel instance
        self.setWindowTitle("Ground Station - Rocket Telemetry")
        self.setGeometry(100, 100, 1400, 900)
        self.setup_ui()
        
        # Start serial connection if port was provided
        if self.selected_port:
            self.start_serial_connection()
    
    def setup_ui(self):
        """Setup the main user interface."""
        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        main_layout = QVBoxLayout(central_widget)
        
        # Control bar at top
        control_layout = QHBoxLayout()
        self.status_label = QLabel("Status: Initializing...")
        control_layout.addWidget(self.status_label)
        control_layout.addStretch()
        
        # Pyro Charges button
        self.pyro_btn = QPushButton("Pyro Charges")
        self.pyro_btn.setStyleSheet("""
            QPushButton {
                background-color: #ff6b35;
                color: #1e1e1e;
                border: none;
                padding: 5px 15px;
                font-weight: bold;
                border-radius: 3px;
            }
            QPushButton:hover {
                background-color: #ff8555;
            }
            QPushButton:pressed {
                background-color: #e55525;
            }
        """)
        self.pyro_btn.clicked.connect(self.open_pyro_panel)
        control_layout.addWidget(self.pyro_btn)
        
        # Add Clear All button
        self.clear_btn = QPushButton("Clear All")
        self.clear_btn.setStyleSheet("""
            QPushButton {
                background-color: #00d4ff;
                color: #1e1e1e;
                border: none;
                padding: 5px 15px;
                font-weight: bold;
                border-radius: 3px;
            }
            QPushButton:hover {
                background-color: #33ddff;
            }
            QPushButton:pressed {
                background-color: #00bbdd;
            }
        """)
        self.clear_btn.clicked.connect(self.clear_all_graphs)
        control_layout.addWidget(self.clear_btn)
        
        self.pause_btn = QPushButton("Pause")
        self.pause_btn.clicked.connect(self.toggle_pause)
        control_layout.addWidget(self.pause_btn)
        
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.clicked.connect(self.stop_connection)
        control_layout.addWidget(self.stop_btn)
        
        main_layout.addLayout(control_layout)
        
        # ALL GRAPHS ON ONE TAB - using grid layout
        graphs_layout = QGridLayout()
        
        # Import and create graph widgets
        self.create_graphs(graphs_layout)
        
        main_layout.addLayout(graphs_layout)
        
        # Status bar at bottom
        self.statusBar().showMessage("Ready")
    
    def create_graphs(self, layout):
        """Create all graphs in a grid layout on one tab."""
        try:
            # Import your graph widgets
            from Frontend.altitude_graph import AltitudeGraph
            from Frontend.temp_graph import TempGraph
            from Frontend.accel_graphlis import AccelGraphLIS
            from Frontend.accel_graphlsm import AccelGraphLSM
            from Frontend.ang_graph import AngGraph
            
            # Create graph instances
            self.altitude_graph = AltitudeGraph()
            self.temp_graph = TempGraph()
            self.accel_lis_graph = AccelGraphLIS()
            self.accel_lsm_graph = AccelGraphLSM()
            self.ang_graph = AngGraph()
            
            # Add to grid layout (2x2 grid)
            # Row 0: Altitude (left), Temperature (right)
            # Row 1: Accel LIS (left), Accel LSM (right)
            layout.addWidget(self.altitude_graph, 0, 0)
            layout.addWidget(self.temp_graph, 0, 1)
            layout.addWidget(self.accel_lis_graph, 1, 0)
            layout.addWidget(self.accel_lsm_graph, 1, 1)
            layout.addWidget(self.ang_graph, 1, 2)
            
        except ImportError as e:
            print(f"Warning: Could not import graph widgets: {e}")
            # Create placeholder labels if graphs not found
            layout.addWidget(QLabel("Altitude Graph - Import Failed"), 0, 0)
            layout.addWidget(QLabel("Temperature Graph - Import Failed"), 0, 1)
            layout.addWidget(QLabel("Accel LIS - Import Failed"), 1, 0)
            layout.addWidget(QLabel("Accel LSM - Import Failed"), 1, 1)
            layout.addWidget(QLabel("Angular Velocity - Import Failed"), 1, 2)
    
    def open_pyro_panel(self):
        """Open the pyro charges control panel."""
        if self.pyro_panel is None:
            from Frontend.pyro_panel import PyroPanel
            self.pyro_panel = PyroPanel(self)
            self.pyro_panel.command_signal.connect(self.send_pyro_command)
        
        self.pyro_panel.show()
        self.pyro_panel.raise_()
        self.pyro_panel.activateWindow()
    
    def send_pyro_command(self, command):
        """Send pyro command via serial."""
        if self.streamer and self.streamer.isRunning():
            self.streamer.write_command(command)
            self.update_status(f"Pyro command sent: {command}")
        else:
            self.update_status("Error: No serial connection active")
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.warning(
                self,
                "Connection Error",
                "Cannot send command: Serial connection is not active"
            )
    
    def clear_all_graphs(self):
        """Clear data from all graphs."""
        if hasattr(self, 'altitude_graph'):
            self.altitude_graph.clear_data()
        if hasattr(self, 'temp_graph'):
            self.temp_graph.clear_data()
        if hasattr(self, 'accel_lis_graph'):
            self.accel_lis_graph.clear_data()
        if hasattr(self, 'accel_lsm_graph'):
            self.accel_lsm_graph.clear_data()
        if hasattr(self, 'ang_graph'):
            self.ang_graph.clear_data()
        self.update_status("All graphs cleared")
    
    def start_serial_connection(self):
        """Start SerialStreamer with the selected port."""
        if not self.selected_port:
            self.update_status("No port selected - running in demo mode")
            return
        
        self.streamer = SerialStreamer(port=self.selected_port)
        self.streamer.new_data.connect(self.handle_new_data)
        self.streamer.status.connect(self.update_status)
        self.streamer.finished.connect(self.on_connection_finished)
        self.streamer.start()
        self.update_status(f"Connecting to {self.selected_port}...")
    
    def handle_new_data(self, data):
        """
        Route incoming data to appropriate graphs.
        Data dictionary contains:
        Time, Temp, Pressure, Alt, Gyro_X, Gyro_Y, Gyro_Z,
        Accel_X1, Accel_Y1, Accel_Z1, Accel_X2, Accel_Y2, Accel_Z2, flight_state
        """
        # Update altitude graph
        if hasattr(self, 'altitude_graph') and data.get('Time') is not None and data.get('Alt') is not None:
            self.altitude_graph.update_data(data.get('Time'), data.get('Alt'))
        
        # Update temperature graph
        if hasattr(self, 'temp_graph') and data.get('Time') is not None and data.get('Temp') is not None:
            self.temp_graph.update_data(data.get('Time'), data.get('Temp'))
        
        # Update LIS accelerometer graph (Accel_X1, Y1, Z1)
        if hasattr(self, 'accel_lis_graph'):
            if all(data.get(k) is not None for k in ['Time', 'Accel_X1', 'Accel_Y1', 'Accel_Z1']):
                self.accel_lis_graph.update_data(
                    data.get('Time'),
                    data.get('Accel_X1'),
                    data.get('Accel_Y1'),
                    data.get('Accel_Z1')
                )
        
        # Update LSM accelerometer graph (Accel_X2, Y2, Z2)
        if hasattr(self, 'accel_lsm_graph'):
            if all(data.get(k) is not None for k in ['Time', 'Accel_X2', 'Accel_Y2', 'Accel_Z2']):
                self.accel_lsm_graph.update_data(
                    data.get('Time'),
                    data.get('Accel_X2'),
                    data.get('Accel_Y2'),
                    data.get('Accel_Z2')
                )
    
        # Update Angular Velocity Graph (Ang_X, Y, Z)
        if hasattr(self, 'ang_graph'):
            if all(data.get(k) is not None for k in ['Time', 'Ang_X2', 'Ang_Y2', 'Ang_Z2']):
                self.ang_graph.update_data(
                    data.get('Time'),
                    data.get('Ang_X2'),
                    data.get('Ang_Y2'),
                    data.get('Ang_Z2')
                )
    
    def update_status(self, message):
        """Update status label and status bar with messages."""
        print(f"[STATUS] {message}")
        self.status_label.setText(f"Status: {message}")
        self.statusBar().showMessage(message)
    
    def toggle_pause(self):
        """Pause or resume data streaming."""
        if not self.streamer:
            return
        
        if self.streamer.paused:
            self.streamer.resume()
            self.pause_btn.setText("Pause")
        else:
            self.streamer.pause()
            self.pause_btn.setText("Resume")
    
    def stop_connection(self):
        """Stop the serial connection."""
        if self.streamer:
            self.streamer.stop()
            self.update_status("Stopping connection...")
    
    def on_connection_finished(self):
        """Handle serial connection closing."""
        self.update_status("Serial connection closed")
        self.pause_btn.setEnabled(False)
        self.stop_btn.setEnabled(False)
    
    def closeEvent(self, event):
        """Handle window closing - stop serial connection."""
        if self.streamer:
            self.streamer.stop()
            self.streamer.wait(1000)  # Wait up to 1 second for thread to finish
        event.accept()
