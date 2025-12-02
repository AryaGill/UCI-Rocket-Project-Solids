from PyQt6.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
                             QFileDialog, QStatusBar, QFrame)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QAction
from Frontend.altitude_graph import AltitudeGraph
from Frontend.temp_graph import TempGraph
from Backend.backend import DataStreamer


class GroundStationWindow(QMainWindow):
    """
    Modern styled main window for rocket telemetry display.
    """
    
    def __init__(self):
        super().__init__()
        self.data_streamer = None
        self.csv_file = None
        
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
        
    def _create_graph_area(self):
        """Create the graph display area with styling."""
        graph_layout = QVBoxLayout()
        graph_layout.setSpacing(15)
        
        # Graph container
        graph_container = QHBoxLayout()
        graph_container.setSpacing(15)
        
        # Altitude graph with frame
        alt_frame = self._create_graph_frame()
        alt_layout = QVBoxLayout(alt_frame)
        self.altitude_graph = AltitudeGraph()
        alt_layout.addWidget(self.altitude_graph)
        graph_container.addWidget(alt_frame)
        
        # Temperature graph with frame
        temp_frame = self._create_graph_frame()
        temp_layout = QVBoxLayout(temp_frame)
        self.temp_graph = TempGraph()
        temp_layout.addWidget(self.temp_graph)
        graph_container.addWidget(temp_frame)
        
        graph_layout.addLayout(graph_container)
        
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
        """Open file dialog to load CSV file."""
        file_name, _ = QFileDialog.getOpenFileName(
            self,
            "Select Telemetry CSV File",
            "",
            "CSV Files (*.csv);;All Files (*)"
        )
        
        if file_name:
            self.csv_file = file_name
            self.start_action.setEnabled(True)
            filename_short = file_name.split('/')[-1]
            self.statusBar.showMessage(f'📄 Loaded: {filename_short}')
            
    def start_stream(self):
        """Start streaming data from CSV."""
        if not self.csv_file:
            return
            
        self.data_streamer = DataStreamer(self.csv_file, delay=0.1)
        self.data_streamer.new_data.connect(self.handle_new_data)
        self.data_streamer.finished.connect(self.stream_finished)
        self.data_streamer.start()
        
        # Update menu actions
        self.start_action.setEnabled(False)
        self.pause_action.setEnabled(True)
        self.stop_action.setEnabled(True)
        
        self.statusBar.showMessage('📡 Streaming telemetry data...')
        
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
        
    def handle_new_data(self, data):
        """Route data to appropriate graph widgets."""
        if 'Time' in data and 'Alt' in data:
            self.altitude_graph.update_data(data['Time'], data['Alt'])
        
        if 'Time' in data and 'Temp' in data:
            self.temp_graph.update_data(data['Time'], data['Temp'])
        
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
