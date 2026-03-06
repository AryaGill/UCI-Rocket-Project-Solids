from PyQt6.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QMessageBox,
                             QGridLayout, QPushButton, QLabel, QSlider, QLineEdit)
from PyQt6.QtCore import (Qt, QTimer, QThread, pyqtSignal, QCoreApplication)
from PyQt6.QtGui import QAction
from Backend.backend import SerialStreamer
from Frontend.flight_state_display import FlightStateDisplay

import csv
import sys
import os
import subprocess

class TTSWorker(QThread):
    """Speaks a message"""
    def __init__(self, message, parent=None):
        super().__init__(parent)
        self.message = message
    
    def run(self):
        try:
            import pyttsx3
            engine = pyttsx3.init()
            engine.setProperty('rate', 160)
            engine.say(self.message)
            engine.runAndWait()
        except Exception as e:
            print(f"[TTS] Error: {e}")

class GroundStationWindow(QMainWindow):
    """
    Main window for the Ground Station application.
    Displays real-time telemetry graphs from the flight computer.
    """
    def __init__(self, port=None, parent=None):
        super().__init__(parent)
        if port is None:
            argv = sys.argv[1:]
            if "--port" in argv:
                idx = argv.index("--port")
                if idx + 1 < len(argv):
                    port = argv[idx + 1]
        self.selected_port = port
        self.streamer = None
        self.pyro_panel = None  # Will hold PyroPanel instance
        self.camera_panel = None  # Will hold CameraPanel instance
        self.max_table = None     # Will hold MaxValuesTable instance
        self.ematch_panel = None
        self.camera_is_on = False   # confirmed state from rocket
        self.camera_pending = None  # "ON" or "OFF" waiting for confirmation
        self.max_points = 100 #Default value, allows us to manually control how many data points we want to see
        self.setWindowTitle("Ground Station - Rocket Telemetry")
        #self.setGeometry(100, 100, 1400, 900)

        self._last_flight_state = None
        self._tts_worker = None

        self.setFixedSize(1500, 1000)
        self.move(100, 100)

        self.setup_ui()

        #CSV attributes specifically for testing, this doesn't affect anything else
        self.csv_rows = []
        self.csv_idx = 0
        self.csv_t0_ms = None
        self.csv_last_ms = None
        self.csv_timer = QTimer(self)
        self.csv_timer.setSingleShot(True)
        self.csv_timer.timeout.connect(self._csv_step)

        self.telemetry_log = []
        
        # Start serial connection if port was provided
        if self.selected_port:
            self.start_serial_connection()
    
    def setup_ui(self):
        # Create menu bar
        menubar = self.menuBar()
        system_menu = menubar.addMenu("System")

        hard_reset_action = QAction("Hard Reset", self)
        hard_reset_action.triggered.connect(self.hard_reset)

        system_menu.addAction(hard_reset_action)

        """Setup the main user interface."""
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        main_layout = QVBoxLayout(central_widget)
        
        # Menu bar
        file_menu = self.menuBar().addMenu("File")

        load_csv_action = QAction("Load CSV...", self)
        load_csv_action.triggered.connect(self.load_csv)
        file_menu.addAction(load_csv_action)

        save_csv_action = QAction("Save CSV...", self)
        save_csv_action.triggered.connect(self.save_csv)
        file_menu.addAction(save_csv_action)

        file_menu.addSeparator()

        clear_action = QAction("Clear All", self)
        clear_action.triggered.connect(self.clear_all_graphs)
        file_menu.addAction(clear_action)

        # Control bar at top
        control_layout = QHBoxLayout()
        self.status_label = QLabel("Status: Initializing...")
        control_layout.addWidget(self.status_label)
        control_layout.addStretch()

        # --- Max points slider ---
        points_label = QLabel("Points: 100")
        control_layout.addWidget(points_label)

        self.points_input = QLineEdit("100")
        self.points_input.setFixedWidth(60)
        self.points_input.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.points_slider = QSlider(Qt.Orientation.Horizontal)
        self.points_slider.setMinimum(10)
        self.points_slider.setMaximum(300)
        self.points_slider.setValue(100)
        self.points_slider.setFixedWidth(150)

        def update_max_points(value):
            self.max_points = value
            points_label.setText(f"Points: {value}")
        
        def slider_changed(value):
            self.max_points = value
            self.points_input.setText(str(value))

        self.points_slider.valueChanged.connect(slider_changed)

        def text_changed():
            text = self.points_input.text()
            if not text.isdigit():
                return

            value = int(text)
            value = max(10, min(1000, value))  # clamp

            self.max_points = value
            self.points_slider.setValue(value)

        self.points_input.editingFinished.connect(text_changed)

        self.points_slider.valueChanged.connect(update_max_points)
        control_layout.addWidget(self.points_input)
        control_layout.addWidget(self.points_slider)
        
        # Add Flight State Display
        try:
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
        self.camera_btn = QPushButton("Camera OFF")
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
        self.camera_btn.clicked.connect(self.toggle_camera)
        control_layout.addWidget(self.camera_btn)

        # Airbrakes Servo Test button
        self.servo_btn = QPushButton("⚙ Airbrakes Test")
        self.servo_btn.setStyleSheet("""
            QPushButton {
                background-color: #a855f7;
                color: #ffffff;
                border: none;
                padding: 5px 15px;
                font-weight: bold;
                border-radius: 3px;
            }
            QPushButton:hover {
                background-color: #bf7fff;
            }
            QPushButton:pressed {
                background-color: #8b3dd4;
            }
        """)
        
        self.servo_btn.clicked.connect(self.test_servo_sequence)
        control_layout.addWidget(self.servo_btn)
        
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

         # ── E-Match status panel ─────────────────────────────────────────────
        try:
            from Frontend.ematch_panel import EMatchPanel
            self.ematch_panel = EMatchPanel()
            main_layout.addWidget(self.ematch_panel)
        except ImportError as e:
            print(f"Warning: Could not import EMatchPanel: {e}")
            self.ematch_panel = None
        
        # ALL GRAPHS ON ONE TAB - using grid layout
        graphs_layout = QGridLayout()
        
        # Import and create graph widgets
        self.create_graphs(graphs_layout)
        
        main_layout.addLayout(graphs_layout)
        
        # Status bar at bottom
        self.statusBar().showMessage("Ready")

        # Announce app startup
        self._tts_worker = TTSWorker("Ground station online")
        self._tts_worker.start()

    def hard_reset(self):
        reply = QMessageBox.warning(
            self, "Hard Reset",
            "This will completely restart the Ground Station.\n\n"
            "All current data and connections will be lost.\n\nContinue?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No
        )
        if reply != QMessageBox.StandardButton.Yes:
            return

        if self.streamer:
            self.streamer.stop()
            self.streamer.wait(2000)

        script_path = self._resolve_entry_script()
        if not script_path:
            QMessageBox.critical(self, "Hard Reset Failed",
                "Could not determine the entry script path.\n\n"
                "Make sure the app is launched as:\n  python main.py")
            return

        cmd = [sys.executable, script_path]
        if self.selected_port:
            cmd += ["--port", self.selected_port]

        try:
            kwargs = {"cwd": os.path.dirname(script_path)}
            if sys.platform == "win32":
                kwargs["creationflags"] = subprocess.DETACHED_PROCESS | subprocess.CREATE_NEW_PROCESS_GROUP
            else:
                kwargs["start_new_session"] = True  # Linux/macOS: detach from parent process group

            subprocess.Popen(cmd, **kwargs)
        except Exception as e:
            QMessageBox.critical(self, "Hard Reset Failed", f"Could not relaunch:\n{e}")
            return

        # Give the new process 500ms to start before this one exits
        QTimer.singleShot(500, QCoreApplication.quit)


    def _resolve_entry_script(self) -> str:
        # 1. Most reliable: set by Python itself when running a .py file directly
        try:
            import __main__
            main_file = getattr(__main__, "__file__", None)
            if main_file:
                path = os.path.abspath(main_file)
                if os.path.isfile(path):
                    return path
        except Exception:
            pass

        # 2. sys.argv[0] as a direct path
        argv = sys.argv[:]
        if argv:
            candidate = os.path.abspath(argv[0])
            if os.path.isfile(candidate):
                return candidate

        # 3. argv tokens split on spaces (paths with spaces)
        for i in range(1, len(argv) + 1):
            candidate = os.path.abspath(" ".join(argv[:i]))
            if os.path.isfile(candidate):
                return candidate

        return ""

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

            # Max values table fills the empty slot: row 0, col 2
            try:
                from Frontend.max_vals import MaxValuesTable
                self.max_table = MaxValuesTable()
                layout.addWidget(self.max_table, 0, 2)
            except ImportError as e:
                print(f"Warning: Could not import MaxValuesTable: {e}")
                self.max_table = None
            
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
    
    '''
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
    '''
    def send_camera_command(self, command):
        """Send camera command via serial and wait for confirmation."""
        if self.streamer and self.streamer.isRunning():
            self.streamer.write_command(command)
            self.camera_pending = command  # store requested state
            self.update_status(f"Camera command sent: {command}")
        else:
            self.update_status("Error: No serial connection active")
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.warning(
                self,
                "Connection Error",
                "Cannot send command: Serial connection is not active"
            )
    
    def toggle_camera(self):
        """Toggle camera on/off via serial command."""
        command = "OFF" if self.camera_is_on else "ON"
        if self.streamer and self.streamer.isRunning():
            self.streamer.write_command(command)
            self.camera_pending = command
            self.update_status(f"Camera command sent: {command}")
        else:
            self.update_status("Error: No serial connection active")
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Connection Error",
                                "Cannot send command: Serial connection is not active")

    def test_servo_sequence(self):
        """Send SERVO SEQUENCE command with confirmation dialog."""
        from PyQt6.QtWidgets import QMessageBox

        reply = QMessageBox.question(
            self,
            "Airbrakes Servo Test",
            "Send SERVO SEQUENCE command?\n\n"
            "The airbrakes servo will run through its full test sequence.\n"
            "Ensure the airbrakes are clear of obstructions before continuing.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No
        )

        if reply != QMessageBox.StandardButton.Yes:
            return

        if self.streamer and self.streamer.isRunning():
            self.streamer.write_command("SERVO SEQUENCE")
            self.update_status("⚙ Airbrakes servo sequence triggered")
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
        if self.max_table is not None:
            self.max_table.reset()
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
        new_state = data.get('flight_state')

        # Handle camera confirmation
        if data.get("camera_status") is not None and self.camera_pending is not None:
            confirmed = int(data.get("camera_status"))

            if confirmed == 1 and self.camera_pending == "ON":
                self.camera_is_on = True
                self.camera_btn.setText("Camera ON")
                self.camera_pending = None
                self.update_status("Camera successfully turned ON")

            elif confirmed == 0 and self.camera_pending == "OFF":
                self.camera_is_on = False
                self.camera_btn.setText("Camera OFF")
                self.camera_pending = None
                self.update_status("Camera successfully turned OFF")

                self.telemetry_log.append(data.copy())
        
        # Update max values table
        if self.max_table is not None:
            self.max_table.update_data(data)

        # Update e-match voltage indicators
        if self.ematch_panel is not None:
            self.ematch_panel.update_data(data)

        # Update flight state display
        if hasattr(self, 'flight_state_display') and data.get('flight_state') is not None:
            new_state = data.get('flight_state')
            self.flight_state_display.update_state(data.get('flight_state'))

            if new_state != self._last_flight_state:
                state_name = FlightStateDisplay.FLIGHT_STATES.get(new_state, "Unknown state")
                print(state_name)
                self._tts_worker = TTSWorker(state_name)
                self._tts_worker.start()
                

            if new_state == 2 and getattr(self, '_last_flight_state', None) != 2:
                self.clear_all_graphs()
                self.update_status("Launch detected - graphs cleared")
            self._last_flight_state = new_state

        # Update altitude graph
        if hasattr(self, 'altitude_graph') and data.get('Time') is not None and data.get('Alt') is not None:
            self.altitude_graph.update_data(data.get('Time'), data.get('Alt'), data.get('Filtered_Alt'), max_points=self.max_points)
        
        # Update temperature graph
        if hasattr(self, 'temp_graph') and data.get('Time') is not None and data.get('Temp') is not None:
            self.temp_graph.update_data(data.get('Time'), data.get('Temp'), max_points=self.max_points)
        
        # Update LIS accelerometer graph (Accel_X1, Y1, Z1)
        if hasattr(self, 'accel_lis_graph'):
            if all(data.get(k) is not None for k in ['Time', 'Accel_X1', 'Accel_Y1', 'Accel_Z1']):
                self.accel_lis_graph.update_data(
                    data.get('Time'),
                    data.get('Accel_X1'),
                    data.get('Accel_Y1'),
                    data.get('Accel_Z1'),
                    max_points=self.max_points
                )

        # Update world accel graph (Accel_world_x, Accel_world_y, Accel_world_z)
        if hasattr(self, 'accel_world_graph'):
            if all(data.get(k) is not None for k in ['Time', 'Accel_world_x', 'Accel_world_y', 'Accel_world_z']):
                self.accel_world_graph.update_data(
                    data.get('Time'),
                    data.get('Accel_world_x'),
                    data.get('Accel_world_y'),
                    data.get('Accel_world_z'),
                    max_points=self.max_points
                )

    
        # Update Angular Velocity Graph (Ang_X, Y, Z)
        if hasattr(self, 'ang_graph'):
            if all(data.get(k) is not None for k in ['Time', 'Gyro_X', 'Gyro_Y', 'Gyro_Z']):
                self.ang_graph.update_data(
                    data.get('Time'),
                    data.get('Gyro_X'),
                    data.get('Gyro_Y'),
                    data.get('Gyro_Z'),
                    max_points=self.max_points
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

    def save_csv(self):
        from PyQt6.QtWidgets import QFileDialog

        if not self.telemetry_log:
            self.update_status("No data to save")
            return

        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Telemetry CSV",
            "telemetry.csv",
            "CSV Files (*.csv)"
        )

        if not path:
            return

        try:
            # Collect all possible keys across packets
            fieldnames = set()
            for row in self.telemetry_log:
                fieldnames.update(row.keys())

            fieldnames = sorted(fieldnames)

            with open(path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(self.telemetry_log)

            self.update_status(f"Saved {len(self.telemetry_log)} rows to {path}")

        except Exception as e:
            self.update_status(f"CSV save error: {e}")

    def load_csv(self):
        from PyQt6.QtWidgets import QFileDialog

        path, _ = QFileDialog.getOpenFileName(
            self, "Open Telemetry CSV", "", "CSV Files (*.csv);;All Files (*)"
        )
        if not path:
            return

        # Stop live streaming
        if self.streamer and self.streamer.isRunning():
            self.streamer.pause()
            self.pause_btn.setText("Resume")

        # Stop any existing CSV playback
        if self.csv_timer.isActive():
            self.csv_timer.stop()

        self.clear_all_graphs()
        self.update_status(f"Loading CSV: {path}")

        def to_float(x):
            try:
                return float(x)
            except Exception:
                return None

        # Load all rows once, then play gradually
        rows = []
        try:
            with open(path, "r", newline="", encoding="utf-8", errors="ignore") as f:
                reader = csv.DictReader(f)
                for r in reader:
                    r = {k.lower(): v for k, v in r.items()}  # ← add this line
                    t_ms = to_float(r.get("time"))
                    if t_ms is None:
                        continue
                    rows.append(r)

            if not rows:
                self.update_status("CSV is empty / no valid rows")
                return

            self.csv_rows = rows
            self.csv_idx = 0
            self.csv_t0_ms = to_float(rows[0].get("time"))
            self.csv_last_ms = self.csv_t0_ms

            self.update_status(f"CSV loaded ({len(rows)} rows). Playing...")
            self._csv_step()  # start playback

        except Exception as e:
            self.update_status(f"CSV load error: {e}")

    def _csv_step(self):
        if self.csv_idx >= len(self.csv_rows):
            self.update_status("CSV playback finished")
            return

        r = self.csv_rows[self.csv_idx]

        def to_float(x):
            try:
                return float(x)
            except Exception:
                return None

        t_ms = to_float(r.get("time"))
        if t_ms is None:
            self.csv_idx += 1
            self.csv_timer.start(0)
            return

        # Convert epoch ms to seconds since start for your graphs
        t_sec = (t_ms - self.csv_t0_ms) / 1000.0

        data = {
            "Time":          t_sec,
            "Alt":           to_float(r.get("alt")),
            "Filtered_Alt":  to_float(r.get("filtered_alt")),   # was None hardcoded

            "Gyro_X":        to_float(r.get("gyro_x")),
            "Gyro_Y":        to_float(r.get("gyro_y")),
            "Gyro_Z":        to_float(r.get("gyro_z")),

            "Accel_X1":      to_float(r.get("accel_x1")),       # was "acc_x"
            "Accel_Y1":      to_float(r.get("accel_y1")),       # was "acc_y"
            "Accel_Z1":      to_float(r.get("accel_z1")),       # was "acc_z"

            "Accel_world_x": to_float(r.get("accel_world_x")),  # was "acc_x_2"
            "Accel_world_y": to_float(r.get("accel_world_y")),  # was "acc_y_2"
            "Accel_world_z": to_float(r.get("accel_world_z")),  # was "acc_z_2"

            "Mag_X":         to_float(r.get("mag_x")),
            "Mag_Y":         to_float(r.get("mag_y")),
            "Mag_Z":         to_float(r.get("mag_z")),

            "Temp":          to_float(r.get("temp")),
            "Pressure":      to_float(r.get("pressure")),
            "flight_state":  to_float(r.get("flight_state")),
        }

        self.handle_new_data(data)
        self.csv_idx += 1

        # Pace next update based on the time difference between consecutive rows
        delay_ms = 0
        if self.csv_idx < len(self.csv_rows):
            next_t_ms = to_float(self.csv_rows[self.csv_idx].get("time"))
            if next_t_ms is not None and self.csv_last_ms is not None:
                delay_ms = int(max(0, min(200, next_t_ms - self.csv_last_ms)))  # cap to keep UI responsive
                self.csv_last_ms = next_t_ms

        self.csv_timer.start(delay_ms)