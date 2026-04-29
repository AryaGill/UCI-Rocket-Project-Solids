from PyQt6.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QMessageBox,
                             QGridLayout, QPushButton, QLabel, QSlider, QLineEdit,
                             QScrollArea,)
from PyQt6.QtCore import (Qt, QTimer, QThread, pyqtSignal, QCoreApplication)
from PyQt6.QtGui import QAction
from Backend.backend import SerialStreamer
from Frontend.flight_state_display import FlightStateDisplay

import csv
import sys
import os
import subprocess
import queue

import queue

class TTSWorker(QThread):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._queue = queue.Queue()
        self._stop = False

    def say(self, message: str):
        self._queue.put(message)

    def stop(self):
        self._stop = True
        self._queue.put(None)

    def run(self):
        while not self._stop:
            try:
                msg = self._queue.get(timeout=1.0)
                if msg is None:
                    break
                if sys.platform == "darwin":
                    subprocess.run(["say", msg])
                elif sys.platform == "win32":
                    subprocess.run(["powershell", "-Command",
                        f'Add-Type -AssemblyName System.Speech; '
                        f'(New-Object System.Speech.Synthesis.SpeechSynthesizer).Speak("{msg}")'])
                else:  # Linux
                    subprocess.run(["espeak", msg])
            except queue.Empty:
                continue
            except Exception as e:
                print(f"[TTS] Error: {e}")

class GroundStationWindow(QMainWindow):
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
        self.pyro_panel = None
        self.camera_panel = None
        self.max_table = None
        self.ematch_panel = None
        self.airbrakes_panel = None
        self.camera_is_on = False
        self.camera_pending = None
        self.max_points = 100
        self.setWindowTitle("Ground Station - Rocket Telemetry")

        self._last_flight_state = None
        self._tts_worker = TTSWorker()
        self._tts_worker.start()

        self.setMinimumSize(1100, 700)
        self.move(100, 100)
        self.showMaximized()

        self.setup_ui()

        self.csv_rows = []
        self.csv_idx = 0
        self.csv_t0_ms = None
        self.csv_last_ms = None
        self.csv_timer = QTimer(self)
        self.csv_timer.setSingleShot(True)
        self.csv_timer.timeout.connect(self._csv_step)

        self.telemetry_log = []
        
        if self.selected_port:
            self.start_serial_connection()

        self._last_announced_m = 0
    
    def setup_ui(self):
        menubar = self.menuBar()
        system_menu = menubar.addMenu("System")

        hard_reset_action = QAction("Hard Reset", self)
        hard_reset_action.triggered.connect(self.hard_reset)
        system_menu.addAction(hard_reset_action)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        main_layout = QVBoxLayout(central_widget)
        
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

        control_layout = QHBoxLayout()
        self.status_label = QLabel("Status: Initializing...")
        control_layout.addWidget(self.status_label)
        control_layout.addStretch()

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
            value = max(10, min(1000, value))
            self.max_points = value
            self.points_slider.setValue(value)

        self.points_input.editingFinished.connect(text_changed)
        self.points_slider.valueChanged.connect(update_max_points)
        control_layout.addWidget(self.points_input)
        control_layout.addWidget(self.points_slider)
        
        try:
            self.flight_state_display = FlightStateDisplay()
            control_layout.addWidget(self.flight_state_display)
        except ImportError as e:
            print(f"Warning: Could not import flight state display: {e}")
        
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
            QPushButton:hover { background-color: #ff8555; }
            QPushButton:pressed { background-color: #e55525; }
        """)
        self.pyro_btn.clicked.connect(self.open_pyro_panel)
        control_layout.addWidget(self.pyro_btn)
        
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
            QPushButton:hover { background-color: #ffff33; }
            QPushButton:pressed { background-color: #e6e200; }
        """)
        self.camera_btn.clicked.connect(self.toggle_camera)
        control_layout.addWidget(self.camera_btn)

        self.gyrocal_btn = QPushButton("🔄 Gyro Cal")
        self.gyrocal_btn.setStyleSheet("""
            QPushButton {
                background-color: #38b0fb;
                color: #1e1e1e;
                border: none;
                padding: 5px 15px;
                font-weight: bold;
                border-radius: 3px;
            }
            QPushButton:hover { background-color: #60c4ff; }
            QPushButton:pressed { background-color: #1a90d4; }
        """)
        self.gyrocal_btn.clicked.connect(self.send_gyrocal)
        control_layout.addWidget(self.gyrocal_btn)

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
            QPushButton:hover { background-color: #bf7fff; }
            QPushButton:pressed { background-color: #8b3dd4; }
        """)
        self.servo_btn.clicked.connect(self.test_servo_sequence)
        control_layout.addWidget(self.servo_btn)
        
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
            QPushButton:hover { background-color: #33ddff; }
            QPushButton:pressed { background-color: #00bbdd; }
        """)
        self.clear_btn.clicked.connect(self.clear_all_graphs)
        control_layout.addWidget(self.clear_btn)
        
        self.pause_btn = QPushButton("Pause")
        self.pause_btn.clicked.connect(self.toggle_pause)
        control_layout.addWidget(self.pause_btn)
        
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.clicked.connect(self.stop_connection)
        control_layout.addWidget(self.stop_btn)

        self.clearSD_btn = QPushButton("Clear SD")
        self.clearSD_btn.clicked.connect(self.send_clearSD)
        control_layout.addWidget(self.clearSD_btn)
        
        main_layout.addLayout(control_layout)

        try:
            from Frontend.ematch_panel import EMatchPanel
            from Frontend.airbrakes_panel import AirbrakesPanel

            self.ematch_panel = EMatchPanel()
            self.airbrakes_panel = AirbrakesPanel()
            top_row = QHBoxLayout()
            top_row.setContentsMargins(0, 0, 0, 0)
            top_row.setSpacing(12)

            top_row.addWidget(self.airbrakes_panel)
            top_row.addStretch()
            top_row.addWidget(self.ematch_panel)

            main_layout.addLayout(top_row)
        except ImportError as e:
            print(f"Warning: Could not import EMatchPanel: {e}")
            self.ematch_panel = None
        
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll_area.setStyleSheet("QScrollArea { border: none; background-color: #1e1e1e; }")

        graphs_container = QWidget()
        graphs_container.setStyleSheet("background-color: #1e1e1e;")
        graphs_layout = QGridLayout(graphs_container)
        graphs_layout.setContentsMargins(4, 4, 4, 4)
        graphs_layout.setSpacing(6)

        self.create_graphs(graphs_layout)

        scroll_area.setWidget(graphs_container)
        main_layout.addWidget(scroll_area)
        
        self.statusBar().showMessage("Ready")
        self._tts_worker.say("Ground station online")

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
                kwargs["start_new_session"] = True

            subprocess.Popen(cmd, **kwargs)
        except Exception as e:
            QMessageBox.critical(self, "Hard Reset Failed", f"Could not relaunch:\n{e}")
            return

        QTimer.singleShot(500, QCoreApplication.quit)

    def _resolve_entry_script(self) -> str:
        try:
            import __main__
            main_file = getattr(__main__, "__file__", None)
            if main_file:
                path = os.path.abspath(main_file)
                if os.path.isfile(path):
                    return path
        except Exception:
            pass

        argv = sys.argv[:]
        if argv:
            candidate = os.path.abspath(argv[0])
            if os.path.isfile(candidate):
                return candidate

        for i in range(1, len(argv) + 1):
            candidate = os.path.abspath(" ".join(argv[:i]))
            if os.path.isfile(candidate):
                return candidate

        return ""

    def create_graphs(self, layout):
        try:
            from Frontend.altitude_graph import AltitudeGraph
            from Frontend.temp_graph import TempGraph
            from Frontend.accel_graphlis import AccelGraphLIS
            from Frontend.accel_graphworld import AccelGraphWorld
            from Frontend.ang_graph import AngGraph
            from Frontend.mag_graph import MagGraph
            from Frontend.rpy_graph import RPYGraph
            from Frontend.velocity_graph import VelocityGraph
            from Frontend.quaternion_graph import QuaternionGraph
            from Frontend.ab_deployment_graph import ABGraph
            from Frontend.apogee_graph import APGraph

            self.altitude_graph = AltitudeGraph()
            self.temp_graph = TempGraph()
            self.accel_lis_graph = AccelGraphLIS()
            self.accel_world_graph = AccelGraphWorld()
            self.ang_graph = AngGraph()
            self.mag_graph = MagGraph()
            self.rpy_graph = RPYGraph()
            self.velocity_graph = VelocityGraph()
            self.quaternion_graph = QuaternionGraph()
            self.ab_graph = ABGraph()
            self.ap_graph = APGraph()

            for graph in (self.altitude_graph, self.temp_graph,
                    self.accel_lis_graph, self.accel_world_graph,
                    self.ang_graph, self.mag_graph,
                    self.rpy_graph, self.velocity_graph,
                    self.quaternion_graph, self.ab_graph,
                    self.ap_graph):
                graph.setMinimumSize(320, 340)

            layout.addWidget(self.altitude_graph,      0, 0)
            layout.addWidget(self.temp_graph,           0, 1)
            layout.addWidget(self.accel_lis_graph,      1, 0)
            layout.addWidget(self.accel_world_graph,    1, 1)
            layout.addWidget(self.ang_graph,            1, 2)
            layout.addWidget(self.mag_graph,            2, 0)
            layout.addWidget(self.rpy_graph,            2, 1)
            layout.addWidget(self.velocity_graph,       2, 2)
            layout.addWidget(self.quaternion_graph,     3, 0)
            layout.addWidget(self.ab_graph,             3, 1)
            layout.addWidget(self.ap_graph,             3, 2)

            for col in range(3):
                layout.setColumnStretch(col, 1)
            for row in range(4):
                layout.setRowStretch(row, 1)

            try:
                from Frontend.max_vals import MaxValuesTable
                self.max_table = MaxValuesTable()
                self.max_table.setMinimumSize(200, 260)
                layout.addWidget(self.max_table, 0, 2)
            except ImportError as e:
                print(f"Warning: Could not import MaxValuesTable: {e}")
                self.max_table = None

        except ImportError as e:
            print(f"Warning: Could not import graph widgets: {e}")
            layout.addWidget(QLabel("Altitude Graph - Import Failed"), 0, 0)
            layout.addWidget(QLabel("Temperature Graph - Import Failed"), 0, 1)
            layout.addWidget(QLabel("Accel LIS - Import Failed"), 1, 0)
            layout.addWidget(QLabel("Accel LSM - Import Failed"), 1, 1)
            layout.addWidget(QLabel("Angular Velocity - Import Failed"), 1, 2)
            layout.addWidget(QLabel("Mag Graph - Import Failed"), 2, 0)
    
    def open_pyro_panel(self):
        if self.pyro_panel is None:
            from Frontend.pyro_panel import PyroPanel
            self.pyro_panel = PyroPanel(self)
            self.pyro_panel.command_signal.connect(self.send_pyro_command)
        self.pyro_panel.show()
        self.pyro_panel.raise_()
        self.pyro_panel.activateWindow()
    
    def open_camera_panel(self):
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
            QMessageBox.warning(self, "Connection Error",
                                "Cannot send command: Serial connection is not active")
    
    def send_camera_command(self, command):
        """Send camera command via serial and wait for confirmation."""
        if self.streamer and self.streamer.isRunning():
            self.streamer.write_command(command)
            self.camera_pending = command
            self.update_status(f"Camera command sent: {command}")
        else:
            self.update_status("Error: No serial connection active")
            from PyQt6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Connection Error",
                                "Cannot send command: Serial connection is not active")
    
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
            self, "Airbrakes Servo Test",
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
            QMessageBox.warning(self, "Connection Error",
                                "Cannot send command: Serial connection is not active")
    
    def send_gyrocal(self):
        """Send GYROCAL command to calibrate the gyroscope."""
        from PyQt6.QtWidgets import QMessageBox
        reply = QMessageBox.question(
            self, "Gyro Calibration",
            "Send GYROCAL command?\n\n"
            "The rocket must be stationary during calibration.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No
        )
        if reply != QMessageBox.StandardButton.Yes:
            return

        if self.streamer and self.streamer.isRunning():
            self.streamer.write_command("GYROCAL")
            self.update_status("🔄 GYROCAL command sent")
        else:
            self.update_status("Error: No serial connection active")
            QMessageBox.warning(self, "Connection Error",
                                "Cannot send command: Serial connection is not active")
    
    def send_clearSD(self):
        if self.streamer and self.streamer.isRunning():
            self.streamer.write_command("RESET_SD")
            self.update_status("🔄 SD Wiped")
        else:
            self.update_status("Error: No serial connection active")
            QMessageBox.warning(self, "Connection Error",
                                "Cannot send command: Serial connection is not active")

    def clear_all_graphs(self):
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
        if hasattr(self, 'mag_graph'):
            self.mag_graph.clear_data()
        if hasattr(self, 'rpy_graph'):
            self.rpy_graph.clear_data()
        if hasattr(self, 'velocity_graph'):
            self.velocity_graph.clear_data()
        if hasattr(self, 'quaternion_graph'):
            self.quaternion_graph.clear_data()
        if self.max_table is not None:
            self.max_table.reset()
        self.update_status("All graphs cleared")
    
    def start_serial_connection(self):
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
        new_state = data.get('flight_state')

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

        if self.max_table is not None:
            self.max_table.update_data(data)

        if self.ematch_panel is not None:
            self.ematch_panel.update_data(data)

        if self.airbrakes_panel is not None:
            self.airbrakes_panel.update_data(data)

        if hasattr(self, 'flight_state_display') and data.get('flight_state') is not None:
            new_state = data.get('flight_state')
            self.flight_state_display.update_state(data.get('flight_state'))

            if new_state != self._last_flight_state:
                state_name = FlightStateDisplay.FLIGHT_STATES.get(new_state, "Unknown state")
                print(state_name + "This should speak")
                self._tts_worker.say(state_name)

            if new_state == 2 and getattr(self, '_last_flight_state', None) != 2:
                self.clear_all_graphs()
                self.update_status("Launch detected - graphs cleared")
            self._last_flight_state = new_state

        if hasattr(self, 'altitude_graph') and data.get('Time') is not None and data.get('Alt') is not None:
            alt_m = data.get('Alt')
            self.altitude_graph.update_data(data.get('Time'), alt_m, data.get('Filtered_Alt'), max_points=self.max_points)
            
            current_m = int(alt_m // 500)
            if alt_m > 0 and current_m != self._last_announced_m:
                self._tts_worker.say(f"{current_m * 500} meters")
                self._last_announced_m = current_m

            '''
            alt_ft = alt_m / 0.3048
            current_kft = int(alt_ft // 1000)
            if current_kft > 0 and current_kft != self._last_announced_kft:
                self._tts_worker.say(f"{current_kft * 1000} feet")
                self._last_announced_kft = current_kft
            '''

        if hasattr(self, 'temp_graph') and data.get('Time') is not None and data.get('Temp') is not None:
            self.temp_graph.update_data(data.get('Time'), data.get('Temp'), max_points=self.max_points)
        
        if hasattr(self, 'accel_lis_graph'):
            if all(data.get(k) is not None for k in ['Time', 'Accel_X1', 'Accel_Y1', 'Accel_Z1']):
                self.accel_lis_graph.update_data(
                    data.get('Time'), data.get('Accel_X1'),
                    data.get('Accel_Y1'), data.get('Accel_Z1'),
                    max_points=self.max_points)

        if hasattr(self, 'accel_world_graph'):
            if all(data.get(k) is not None for k in ['Time', 'Accel_world_x', 'Accel_world_y', 'Accel_world_z']):
                self.accel_world_graph.update_data(
                    data.get('Time'), data.get('Accel_world_x'),
                    data.get('Accel_world_y'), data.get('Accel_world_z'),
                    max_points=self.max_points)

        if hasattr(self, 'ang_graph'):
            if all(data.get(k) is not None for k in ['Time', 'Gyro_X', 'Gyro_Y', 'Gyro_Z']):
                self.ang_graph.update_data(
                    data.get('Time'), data.get('Gyro_X'),
                    data.get('Gyro_Y'), data.get('Gyro_Z'),
                    max_points=self.max_points)

        if hasattr(self, 'mag_graph'):
            if all(data.get(k) is not None for k in ['Time', 'mag_r', 'mag_p', 'mag_y']):
                self.mag_graph.update_data(
                    data['Time'], data['mag_r'], data['mag_p'], data['mag_y'])
        
        if hasattr(self, 'rpy_graph'):
            if all(data.get(k) is not None for k in ['Time', 'roll', 'pitch', 'yaw']):
                self.rpy_graph.update_data(
                    data['Time'], data['roll'], data['pitch'], data['yaw'],
                    max_points=self.max_points)

        if hasattr(self, 'velocity_graph'):
            if all(data.get(k) is not None for k in ['Time', 'velocity_x', 'velocity_y', 'velocity_z']):
                self.velocity_graph.update_data(
                    data['Time'], data['velocity_x'], data['velocity_y'], data['velocity_z'],
                    max_points=self.max_points)

        if hasattr(self, 'quaternion_graph'):
            if all(data.get(k) is not None for k in ['Time', 'Quaternion_W', 'Quaternion_X', 'Quaternion_Y', 'Quaternion_Z']):
                self.quaternion_graph.update_data(
                    data['Time'], data['Quaternion_W'], data['Quaternion_X'],
                    data['Quaternion_Y'], data['Quaternion_Z'],
                    max_points=self.max_points)
        
        if hasattr(self, 'ab_graph'):
            if all(data.get(k) is not None for k in ['Time', 'AB_Deployment']):
                self.ab_graph.update_data(
                        data['Time'], (data['AB_Deployment'] / 63) * 100, max_points = self.max_points)
                
        if hasattr(self, 'ap_graph'):
            if all(data.get(k) is not None for k in ['Time', 'pred_apo']):
                self.ap_graph.update_data(
                        data['Time'], data['pred_apo'], max_points = self.max_points)

    
    def update_status(self, message):
        print(f"[STATUS] {message}")
        self.status_label.setText(f"Status: {message}")
        self.statusBar().showMessage(message)
    
    def toggle_pause(self):
        if not self.streamer:
            return
        if self.streamer.paused:
            self.streamer.resume()
            self.pause_btn.setText("Pause")
        else:
            self.streamer.pause()
            self.pause_btn.setText("Resume")
    
    def stop_connection(self):
        if self.streamer:
            self.streamer.stop()
            self.update_status("Stopping connection...")
    
    def on_connection_finished(self):
        self.update_status("Serial connection closed")
        self.pause_btn.setEnabled(False)
        self.stop_btn.setEnabled(False)
    
    def closeEvent(self, event):
        if self.streamer:
            self.streamer.stop()
            self.streamer.wait(1000)
        self._tts_worker.stop()
        self._tts_worker.wait(2000)
        event.accept()

    def arm_rocket(self):
        from PyQt6.QtWidgets import QMessageBox
        
        reply1 = QMessageBox.warning(
            self, "⚠️ ARM ROCKET - FIRST CONFIRMATION",
            "You are about to ARM the rocket.\n\n"
            "This will enable pyrotechnic charges and prepare the flight computer for launch.\n\n"
            "Are you sure you want to continue?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No
        )
        if reply1 != QMessageBox.StandardButton.Yes:
            return
        
        reply2 = QMessageBox.critical(
            self, "🚀 ARM ROCKET - FINAL CONFIRMATION",
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
                print(f"_ser={self.streamer._ser}, is_open={getattr(self.streamer._ser, 'is_open', 'N/A')}")
                self.update_status("🚀 ARM command sent - ROCKET ARMED")
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
                self.arm_btn.setEnabled(False)
            else:
                self.update_status("Error: No serial connection active")
                QMessageBox.warning(self, "Connection Error",
                                    "Cannot send ARM command: Serial connection is not active")

    def save_csv(self):
        from PyQt6.QtWidgets import QFileDialog

        if not self.telemetry_log:
            self.update_status("No data to save")
            return

        path, _ = QFileDialog.getSaveFileName(
            self, "Save Telemetry CSV", "telemetry.csv", "CSV Files (*.csv)")
        if not path:
            return

        try:
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
            self, "Open Telemetry CSV", "", "CSV Files (*.csv);;All Files (*)")
        if not path:
            return

        if self.streamer and self.streamer.isRunning():
            self.streamer.pause()
            self.pause_btn.setText("Resume")

        if self.csv_timer.isActive():
            self.csv_timer.stop()

        self.clear_all_graphs()
        self.update_status(f"Loading CSV: {path}")

        def to_float(x):
            try:
                return float(x)
            except Exception:
                return None

        rows = []
        try:
            with open(path, "r", newline="", encoding="utf-8", errors="ignore") as f:
                reader = csv.DictReader(f)
                for r in reader:
                    r = {k.lower(): v for k, v in r.items()}
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
            self._csv_step()

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

        t_sec = (t_ms - self.csv_t0_ms) / 1000.0

        data = {
            "Time":          t_sec,
            "Alt":           to_float(r.get("alt")),
            "Filtered_Alt":  to_float(r.get("filtered_alt")),
            "Gyro_X":        to_float(r.get("gyro_x")),
            "Gyro_Y":        to_float(r.get("gyro_y")),
            "Gyro_Z":        to_float(r.get("gyro_z")),
            "Accel_X1":      to_float(r.get("accel_x1")),
            "Accel_Y1":      to_float(r.get("accel_y1")),
            "Accel_Z1":      to_float(r.get("accel_z1")),
            "Accel_world_x": to_float(r.get("accel_world_x")),
            "Accel_world_y": to_float(r.get("accel_world_y")),
            "Accel_world_z": to_float(r.get("accel_world_z")),
            "mag_r":         to_float(r.get("mag_r")),
            "mag_p":         to_float(r.get("mag_p")),
            "mag_y":         to_float(r.get("mag_y")),
            "Temp":          to_float(r.get("temp")),
            "Pressure":      to_float(r.get("pressure")),
            "flight_state":  to_float(r.get("flight_state")),
        }

        self.handle_new_data(data)
        self.csv_idx += 1

        delay_ms = 0
        if self.csv_idx < len(self.csv_rows):
            next_t_ms = to_float(self.csv_rows[self.csv_idx].get("time"))
            if next_t_ms is not None and self.csv_last_ms is not None:
                delay_ms = int(max(0, min(200, next_t_ms - self.csv_last_ms)))
                self.csv_last_ms = next_t_ms

        self.csv_timer.start(delay_ms)