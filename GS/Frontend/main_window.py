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
        self.camera_panel = None  # Will hold CameraPanel instance
        self.setWindowTitle("Ground Station - Rocket Telemetry")
        self.setGeometry(100, 100, 1400, 900)
        self.setup_ui()
        
        # Start serial connection if port was provided
        if self.selected_port:
            self.start_serial_connection()
    
    def setup_ui(self):
        """Setup the main user interface."""
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        main_layout = QVBoxLayout(central_widget)
        
        # Control bar at top
        control_layout = QHBoxLayout()
        self.status_label = QLabel("Status: Initializing...")
        control_layout.addWidget(self.status_label)
        control_layout.addStretch()
        
        # Add Flight State Display
        try:
            from Frontend.flight_state_display import FlightStateDisplay
            self.flight_state_display = FlightStateDisplay()
            control_layout.addWidget(self.flight_state_display)
        except ImportError as e:
            print(f"Warning: Could not import flight state display: {e}")
        
        # ARM button with safety lock
        self.arm_btn = QPushButton("🔒 ARM ROCKET")
        self.arm_btn.setMinimumHeight(40)
        self.arm_btn.setStyleSheet("""
            QPushButton {
                background-color: #ff0000;
                color: #ffffff;
                border: 3px solid #ff0000;
                padding: 8px 20px;
                font-weight: bold;
                font-size: 14px;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #ff3333;
                border: 3px solid #ff3333;
            }
            QPushButton:pressed {
                background-color: #cc0000;
            }
        """)
        self.arm_btn.clicked.connect(self.arm_rocket)
        control_layout.addWidget(self.arm_btn)
        
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
        
        # Camera button
        self.camera_btn = QPushButton("Camera")
        self.camera_btn.setStyleSheet("""
            QPushButton {
                background-color: #fffb00;
                color: #1e1e1e;
                border: none;
                padding: 5px 15px;
                font-weight: bold;
                border-radius: 3px;
            }
            QPushButton:hover {
                background-color: #ffff33;
            }
            QPushButton:pressed {
                background-color: #e6e200;
            }
        """)
        self.camera_btn.clicked.connect(self.open_camera_panel)
        control_layout.addWidget(self.camera_btn)
        
        # Clear All button
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
            from Frontend.accel_graphworld import AccelGraphWorld
            from Frontend.ang_graph import AngGraph
            
            # Create graph instances
            self.altitude_graph = AltitudeGraph()
            self.temp_graph = TempGraph()
            self.accel_lis_graph = AccelGraphLIS()
            self.accel_world_graph = AccelGraphWorld()
            self.ang_graph = AngGraph()
            
            # Add to grid layout (2x2 grid)
            # Row 0: Altitude (left), Temperature (right)
            # Row 1: Accel LIS (left), Accel LSM (right)
            layout.addWidget(self.altitude_graph, 0, 0)
            layout.addWidget(self.temp_graph, 0, 1)
            layout.addWidget(self.accel_lis_graph, 1, 0)
            layout.addWidget(self.accel_world_graph, 1, 1)
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
    
    def open_camera_panel(self):
        """Open the camera control panel."""
        if self.camera_panel is None:
            from Frontend.camera_panel import CameraPanel
            self.camera_panel = CameraPanel(self)
            self.camera_panel.command_signal.connect(self.send_camera_command)
        
        self.camera_panel.show()
        self.camera_panel.raise_()
        self.camera_panel.activateWindow()
    
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
    
    def send_camera_command(self, command):
        """Send camera command via serial."""
        if self.streamer and self.streamer.isRunning():
            self.streamer.write_command(command)
            self.update_status(f"Camera command sent: {command}")
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
        if hasattr(self, 'accel_world_graph'):
            self.accel_world_graph.clear_data()
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
        # Update flight state display
        if hasattr(self, 'flight_state_display') and data.get('flight_state') is not None:
            self.flight_state_display.update_state(data.get('flight_state'))
        
        # Update altitude graph
        if hasattr(self, 'altitude_graph') and data.get('Time') is not None and data.get('Alt') is not None:
            self.altitude_graph.update_data(data.get('Time'), data.get('Alt'), data.get('Filtered_Alt'))
        
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

        # Update world accel graph (Accel_world_x, Accel_world_y, Accel_world_z)
        if hasattr(self, 'accel_world_graph'):
            if all(data.get(k) is not None for k in ['Time', 'Accel_world_x', 'Accel_world_y', 'Accel_world_z']):
                self.accel_world_graph.update_data(
                    data.get('Time'),
                    data.get('Accel_world_x'),
                    data.get('Accel_world_y'),
                    data.get('Accel_world_z')
                )

    
        # Update Angular Velocity Graph (Ang_X, Y, Z)
        if hasattr(self, 'ang_graph'):
            if all(data.get(k) is not None for k in ['Time', 'Gyro_X', 'Gyro_Y', 'Gyro_Z']):
                self.ang_graph.update_data(
                    data.get('Time'),
                    data.get('Gyro_X'),
                    data.get('Gyro_Y'),
                    data.get('Gyro_Z')
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

    def arm_rocket(self):
        """Send ARM command with safety confirmation."""
        from PyQt6.QtWidgets import QMessageBox
        
        # First confirmation
        reply1 = QMessageBox.warning(
            self,
            "⚠️ ARM ROCKET - FIRST CONFIRMATION",
            "You are about to ARM the rocket.\n\n"
            "This will enable pyrotechnic charges and prepare the flight computer for launch.\n\n"
            "Are you sure you want to continue?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No
        )
        
        if reply1 != QMessageBox.StandardButton.Yes:
            return
        
        # Second confirmation
        reply2 = QMessageBox.critical(
            self,
            "🚀 ARM ROCKET - FINAL CONFIRMATION",
            "FINAL WARNING!\n\n"
            "Arming the rocket will:\n"
            "• Enable all pyrotechnic circuits\n"
            "• Activate flight detection algorithms\n"
            "• Begin autonomous flight operations\n\n"
            "Confirm ARM command?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No
        )
        
        if reply2 == QMessageBox.StandardButton.Yes:
            if self.streamer and self.streamer.isRunning():
                self.streamer.write_command("ARM")
                self.update_status("🚀 ARM command sent - ROCKET ARMED")
                
                # Change button appearance after arming
                self.arm_btn.setText("✓ ARMED")
                self.arm_btn.setStyleSheet("""
                    QPushButton {
                        background-color: #28dc5e;
                        color: #1e1e1e;
                        border: 3px solid #28dc5e;
                        padding: 8px 20px;
                        font-weight: bold;
                        font-size: 14px;
                        border-radius: 5px;
                    }
                """)
                self.arm_btn.setEnabled(False)  # Disable after arming
            else:
                self.update_status("Error: No serial connection active")
                QMessageBox.warning(
                    self,
                    "Connection Error",
                    "Cannot send ARM command: Serial connection is not active"
                )
