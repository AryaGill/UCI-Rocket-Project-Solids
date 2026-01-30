from PyQt6.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
                             QFileDialog, QStatusBar, QFrame)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QAction
from Frontend.altitude_graph import AltitudeGraph
from Frontend.temp_graph import TempGraph
from Frontend.mag_graph import MagGraph
from Frontend.ang_graph import AngGraph
from Frontend.accel_graphlis import AccelGraphLIS
from Frontend.accel_graphlsm import AccelGraphLSM

from Frontend.pyrodrogue_panel import PyroDroguePanel
from Frontend.status import StatusIndicator

from Backend.backend import DataStreamer, SerialStreamer


class GroundStationWindow(QMainWindow):
    """
    Modern styled main window for rocket telemetry display.
    """
    
    def __init__(self):
        super().__init__()
        self.data_streamer = None
        self.selected_port = None
        self.csv_file = None

        self._t0_ms = None
        
        self.init_ui()
        
    def init_ui(self):
        """Initialize the user interface with modern styling."""
        self.setWindowTitle('🚀 Rocket Ground Station - Telemetry Display')
        self.setGeometry(100, 100, 1500, 900)
        
        # Create menu bar
        self._create_menu_bar()
        
        # Create central widget and main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        main_layout.setSpacing(15)
        main_layout.setContentsMargins(20, 20, 20, 20)
        
        # Graph area (no header)
        graph_layout = self._create_graph_area()
        main_layout.addLayout(graph_layout)

        # Pyro Drogue panel (hidden by default)
        self.pyro_panel = PyroDroguePanel()

        self.pyro_panel.arm_clicked.connect(lambda: self.send_rf_command("MAIN_PRIMARY"))
        self.pyro_panel.fire_clicked.connect(lambda: self.send_rf_command("MAIN_SECONDARY"))
        self.pyro_panel.disarm_clicked.connect(lambda: self.send_rf_command("DROGUE_PRIMARY"))
        self.pyro_panel.abort_clicked.connect(lambda: self.send_rf_command("DROGUE_SECONDARY"))


        main_layout.addWidget(self.pyro_panel)
        self.pyro_panel.setMaximumWidth(self.width() // 2)

        # Right: status indicator
        self.status_indicator = StatusIndicator(diameter=30)
        main_layout.addWidget(self.status_indicator)

        # Example: update status (in-flight)
        self.status_indicator.set_status(True)
        
        # Styled status bar
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        self.statusBar.setStyleSheet("""
            QStatusBar {
                background-color: #2d2d2d;
                color: #00d4ff;
                font-size: 12px;
                padding: 5px;
                border-top: 2px solid #404040;
            }
        """)
        self.statusBar.showMessage('⚡ Ready - Load CSV file to start streaming')
        
    def _create_menu_bar(self):
        """Create menu bar with File and Control menus."""
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu('File')
        
        # Load action
        load_action = QAction('📁 Load CSV...', self)
        load_action.setShortcut('Ctrl+O')
        load_action.triggered.connect(self.load_csv)
        file_menu.addAction(load_action)
        
        file_menu.addSeparator()
        
        # Clear action
        clear_action = QAction('🗑 Clear Graphs', self)
        clear_action.setShortcut('Ctrl+K')
        clear_action.triggered.connect(self.clear_graphs)
        file_menu.addAction(clear_action)
        
        file_menu.addSeparator()
        
        # Exit action
        exit_action = QAction('Exit', self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # Control menu
        control_menu = menubar.addMenu('Control')
        
        # Start action
        self.start_action = QAction('▶ Start Streaming', self)
        self.start_action.setShortcut('Ctrl+S')
        self.start_action.setEnabled(False)
        self.start_action.triggered.connect(self.start_stream)
        control_menu.addAction(self.start_action)
        
        # Pause action
        self.pause_action = QAction('⏸ Pause', self)
        self.pause_action.setShortcut('Ctrl+P')
        self.pause_action.setEnabled(False)
        self.pause_action.triggered.connect(self.toggle_pause)
        control_menu.addAction(self.pause_action)
        
        # Stop action
        self.stop_action = QAction('⏹ Stop Streaming', self)
        self.stop_action.setShortcut('Ctrl+T')
        self.stop_action.setEnabled(False)
        self.stop_action.triggered.connect(self.stop_stream)
        control_menu.addAction(self.stop_action)

        pyro_menu = menubar.addMenu('Pyro Drogue')

        self.toggle_pyro_action = QAction('Show Pyro Controls', self)
        self.toggle_pyro_action.setCheckable(True)
        self.toggle_pyro_action.triggered.connect(self.toggle_pyro_panel)

        pyro_menu.addAction(self.toggle_pyro_action)

    def toggle_pyro_panel(self, checked):
        self.pyro_panel.setVisible(checked)

        if checked:
            self.toggle_pyro_action.setText('Hide Pyro Controls')
            self.statusBar.showMessage('🔥 Pyro Drogue controls enabled')
        else:
            self.toggle_pyro_action.setText('Show Pyro Controls')

    def _create_graph_area(self):
        """Create the graph display area with styling."""
        graph_layout = QVBoxLayout()
        graph_layout.setSpacing(15)
        
        # Graph container
        row1 = QHBoxLayout()
        row2 = QHBoxLayout()
        row1.setSpacing(15)
        row2.setSpacing(15)

        
        # Altitude graph with frame
        alt_frame = self._create_graph_frame()
        alt_layout = QVBoxLayout(alt_frame)
        self.altitude_graph = AltitudeGraph()
        alt_layout.addWidget(self.altitude_graph)
        row1.addWidget(alt_frame)
        
        # Temperature graph with frame
        temp_frame = self._create_graph_frame()
        temp_layout = QVBoxLayout(temp_frame)
        self.temp_graph = TempGraph()
        temp_layout.addWidget(self.temp_graph)
        row1.addWidget(temp_frame)

        # Magnetic Field graph with frame
        mag_frame = self._create_graph_frame()
        mag_layout = QVBoxLayout(mag_frame)
        self.mag_graph = MagGraph()
        mag_layout.addWidget(self.mag_graph)
        row1.addWidget(mag_frame)

        # Angular Speed graph with frame
        ang_frame = self._create_graph_frame()
        ang_layout = QVBoxLayout(ang_frame)
        self.ang_graph = AngGraph()
        ang_layout.addWidget(self.ang_graph)
        row2.addWidget(ang_frame)

        # Angular Speed graph with frame
        accelLIS_frame = self._create_graph_frame()
        accelLIS_layout = QVBoxLayout(accelLIS_frame)
        self.accelLIS_graph = AccelGraphLIS()
        accelLIS_layout.addWidget(self.accelLIS_graph)
        row2.addWidget(accelLIS_frame)

        # Angular Speed graph with frame
        accelLSM_frame = self._create_graph_frame()
        accelLSM_layout = QVBoxLayout(accelLSM_frame)
        self.accelLSM_graph = AccelGraphLSM()
        accelLSM_layout.addWidget(self.accelLSM_graph)
        row2.addWidget(accelLSM_frame)
        
        graph_layout.addLayout(row1)
        graph_layout.addLayout(row2)
        
        return graph_layout
        
    def _create_graph_frame(self):
        """Create a styled frame for graphs."""
        frame = QFrame()
        frame.setStyleSheet("""
            QFrame {
                background-color: #2d2d2d;
                border-radius: 10px;
                border: 2px solid #404040;
            }
        """)
        return frame
        
    def load_csv(self):
        file_name, _ = QFileDialog.getOpenFileName(
            self, "Select Telemetry CSV File", "", "CSV Files (*.csv);;All Files (*)"
        )

        if file_name:
            self.csv_file = file_name
            self.selected_port = None  # <-- ADD THIS LINE (force CSV playback)
            self.start_action.setEnabled(True)
            filename_short = file_name.split('/')[-1]
            self.statusBar.showMessage(f'📄 Loaded: {filename_short}')
            self.start_stream()


    #This defines what the pyro drogue buttons actually do
    def send_rf_command(self, cmd: str):
        """
        if not self.data_streamer:
            self.statusBar.showMessage("No active link.")
            return

        # Works only for SerialStreamer; CSV streamer can’t transmit.
        if hasattr(self.data_streamer, "write_command"):
            self.data_streamer.write_command(cmd)
            self.statusBar.showMessage(f"Sent: {cmd}")
        else:
            self.statusBar.showMessage("Streaming from CSV; cannot transmit.")
        """
        print(f"[UI COMMAND] {cmd}")
        self.statusBar.showMessage(f"Console: {cmd}")

    def start_stream(self):
        """Start streaming data from selected source (serial preferred, CSV fallback)."""
        
        # Prefer serial if selected
        if self.selected_port:
            self.data_streamer = SerialStreamer(self.selected_port, baud=115200)
            self.data_streamer.new_data.connect(self.handle_new_data)
            self.data_streamer.finished.connect(self.stream_finished)
            self.data_streamer.status.connect(self.statusBar.showMessage)
            self.data_streamer.start()
            self.statusBar.showMessage(f"📡 Streaming from {self.selected_port}...")
        else:
            # fallback to CSV (your existing behavior)
            if not self.csv_file:
                return
            self.data_streamer = DataStreamer(self.csv_file, delay=0.1)
            self.data_streamer.new_data.connect(self.handle_new_data)
            self.data_streamer.finished.connect(self.stream_finished)
            self.data_streamer.start()
            self.statusBar.showMessage("📡 Streaming telemetry data (CSV)...")
        self.start_action.setEnabled(False)
        self.pause_action.setEnabled(True)
        self.stop_action.setEnabled(True)

        
    def toggle_pause(self):
        """Pause or resume the data stream."""
        if not self.data_streamer:
            return
            
        if self.pause_action.text() == '⏸ Pause':
            self.data_streamer.pause()
            self.pause_action.setText('▶ Resume')
            self.statusBar.showMessage('⏸ Paused')
        else:
            self.data_streamer.resume()
            self.pause_action.setText('⏸ Pause')
            self.statusBar.showMessage('📡 Streaming telemetry data...')
            
    def stop_stream(self):
        """Stop the data stream."""
        if self.data_streamer:
            self.data_streamer.stop()
            self.data_streamer.wait()
        self.stream_finished()
        
    def stream_finished(self):
        """Handle stream completion."""
        self.start_action.setEnabled(True)
        self.pause_action.setEnabled(False)
        self.pause_action.setText('⏸ Pause')
        self.stop_action.setEnabled(False)
        self.statusBar.showMessage('✅ Stream finished')
        
    def handle_new_data(self, data: dict):
        # normalize keys (lowercase) so CSV + serial can both work
        d = {str(k).strip().lower(): v for k, v in data.items()}

        # time in your CSV is epoch milliseconds; convert to seconds since start
        if "time" not in d:
            return

        t_ms = float(d["time"])
        if not hasattr(self, "_t0_ms") or self._t0_ms is None:
            self._t0_ms = t_ms
        t = (t_ms - self._t0_ms) / 1000.0  # seconds since start

        # altitude
        if "alt" in d:
            self.altitude_graph.update_data(t, float(d["alt"]))

        # magnetometer
        if all(k in d for k in ("mag_x", "mag_y", "mag_z")):
            self.mag_graph.update_data(t, float(d["mag_x"]), float(d["mag_y"]), float(d["mag_z"]))

        # gyro (angular rate)
        if all(k in d for k in ("gyro_x", "gyro_y", "gyro_z")):
            self.ang_graph.update_data(t, float(d["gyro_x"]), float(d["gyro_y"]), float(d["gyro_z"]))

        # accel set 1
        if all(k in d for k in ("acc_x", "acc_y", "acc_z")):
            self.accelLIS_graph.update_data(t, float(d["acc_x"]), float(d["acc_y"]), float(d["acc_z"]))

        # accel set 2
        if all(k in d for k in ("acc_x_2", "acc_y_2", "acc_z_2")):
            self.accelLSM_graph.update_data(t, float(d["acc_x_2"]), float(d["acc_y_2"]), float(d["acc_z_2"]))
        


        
    def clear_graphs(self):
        """Clear all graph data."""
        self.altitude_graph.clear_data()
        self.temp_graph.clear_data()
        self.statusBar.showMessage('🗑 Graphs cleared')
        
    def closeEvent(self, event):
        """Handle application close event."""
        if self.data_streamer and self.data_streamer.isRunning():
            self.data_streamer.stop()
            self.data_streamer.wait()
        event.accept()
