"""
Main GUI for real-time sensor data visualization
Uses PyQt6 and pyqtgraph for high-performance plotting
"""

import sys
from collections import deque
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QPushButton, QLabel, QComboBox, 
                             QMessageBox, QGroupBox)
from PyQt6.QtCore import QTimer, Qt
import pyqtgraph as pg
from serial_backend import SerialBackend
from sensor_config import SENSORS, AXIS_COLORS
import serial.tools.list_ports


class SensorGraph(QWidget):
    """Reusable graph widget for sensor data"""

    def __init__(self, title, columns, y_label, colors, window_size=100):
        super().__init__()
        self.columns = columns
        self.window_size = window_size
        self.colors = colors

        # Data buffers (rolling window)
        self.data_buffers = {col: deque(maxlen=window_size) for col in columns}
        self.time_buffer = deque(maxlen=window_size)
        self.time_counter = 0

        # Setup layout
        layout = QVBoxLayout()
        self.setLayout(layout)

        # Create plot widget
        self.plot_widget = pg.PlotWidget(title=title)
        self.plot_widget.setLabel('left', y_label)
        self.plot_widget.setLabel('bottom', 'Sample')
        self.plot_widget.showGrid(x=True, y=True, alpha=0.3)
        self.plot_widget.addLegend()

        # Create plot curves for each column
        self.curves = {}
        for col, color in zip(columns, colors):
            pen = pg.mkPen(color=color, width=2)
            self.curves[col] = self.plot_widget.plot(
                [], [], 
                pen=pen, 
                name=col.capitalize()
            )

        layout.addWidget(self.plot_widget)

    def update_data(self, data_dict):
        """Update graph with new data point"""
        self.time_buffer.append(self.time_counter)
        self.time_counter += 1

        for col in self.columns:
            if col in data_dict:
                self.data_buffers[col].append(data_dict[col])

        # Update all curves
        time_array = list(self.time_buffer)
        for col, curve in self.curves.items():
            data_array = list(self.data_buffers[col])
            if len(data_array) == len(time_array):
                curve.setData(time_array, data_array)

    def clear_data(self):
        """Clear all data buffers"""
        for buffer in self.data_buffers.values():
            buffer.clear()
        self.time_buffer.clear()
        self.time_counter = 0

        for curve in self.curves.values():
            curve.setData([], [])


class SensorDashboard(QMainWindow):
    """Main dashboard window"""

    def __init__(self, sensor_type='LPS22HHTR'):
        super().__init__()
        self.sensor_type = sensor_type
        self.sensor_config = SENSORS[sensor_type]
        self.backend = SerialBackend(sensor_type=sensor_type)
        self.is_streaming = False

        self.init_ui()

        # Setup update timer
        self.update_timer = QTimer()
        self.update_timer.timeout.connect(self.update_graphs)

    def init_ui(self):
        """Initialize the user interface"""
        self.setWindowTitle(f'Sensor Dashboard - {self.sensor_type}')
        self.setGeometry(100, 100, 1400, 800)

        # Main widget and layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QVBoxLayout()
        main_widget.setLayout(main_layout)

        # Control panel
        control_group = self.create_control_panel()
        main_layout.addWidget(control_group)

        # Graph container
        graph_layout = QVBoxLayout()

        # Create graphs based on sensor configuration
        self.graphs = []
        for graph_config in self.sensor_config['graphs']:
            graph = SensorGraph(
                title=graph_config['title'],
                columns=graph_config['columns'],
                y_label=graph_config['y_label'],
                colors=graph_config['colors'],
                window_size=self.sensor_config['window_size']
            )
            self.graphs.append(graph)
            graph_layout.addWidget(graph)

        main_layout.addLayout(graph_layout)

        # Status bar
        self.statusBar().showMessage('Ready - Not connected')

    def create_control_panel(self):
        """Create the control panel with buttons and port selection"""
        control_group = QGroupBox("Control Panel")
        control_layout = QHBoxLayout()

        # Port selection
        port_label = QLabel("Serial Port:")
        control_layout.addWidget(port_label)

        self.port_combo = QComboBox()
        self.refresh_ports()
        control_layout.addWidget(self.port_combo)

        # Refresh button
        self.refresh_btn = QPushButton("Refresh Ports")
        self.refresh_btn.clicked.connect(self.refresh_ports)
        control_layout.addWidget(self.refresh_btn)

        # Connect button
        self.connect_btn = QPushButton("Connect")
        self.connect_btn.clicked.connect(self.toggle_connection)
        control_layout.addWidget(self.connect_btn)

        # Start/Stop streaming button
        self.stream_btn = QPushButton("Start Streaming")
        self.stream_btn.clicked.connect(self.toggle_streaming)
        self.stream_btn.setEnabled(False)
        control_layout.addWidget(self.stream_btn)

        # Clear button
        self.clear_btn = QPushButton("Clear Graphs")
        self.clear_btn.clicked.connect(self.clear_graphs)
        control_layout.addWidget(self.clear_btn)

        control_layout.addStretch()
        control_group.setLayout(control_layout)

        return control_group

    def refresh_ports(self):
        """Refresh the list of available serial ports"""
        self.port_combo.clear()
        ports = serial.tools.list_ports.comports()

        for port in ports:
            self.port_combo.addItem(f"{port.device} - {port.description}", port.device)

        if self.port_combo.count() == 0:
            self.port_combo.addItem("No ports found", None)

    def toggle_connection(self):
        """Connect or disconnect from serial port"""
        if not self.backend.serial_conn or not self.backend.serial_conn.is_open:
            # Connect
            port = self.port_combo.currentData()

            if port is None:
                QMessageBox.warning(self, "Error", "No serial port selected")
                return

            if self.backend.connect(port):
                self.connect_btn.setText("Disconnect")
                self.stream_btn.setEnabled(True)
                self.port_combo.setEnabled(False)
                self.refresh_btn.setEnabled(False)
                self.statusBar().showMessage(f'Connected to {port}')
            else:
                QMessageBox.critical(self, "Connection Error", f"Failed to connect to {port}")
        else:
            # Disconnect
            if self.is_streaming:
                self.toggle_streaming()

            self.backend.disconnect()
            self.connect_btn.setText("Connect")
            self.stream_btn.setEnabled(False)
            self.port_combo.setEnabled(True)
            self.refresh_btn.setEnabled(True)
            self.statusBar().showMessage('Disconnected')

    def toggle_streaming(self):
        """Start or stop data streaming"""
        if not self.is_streaming:
            # Start streaming
            self.is_streaming = True
            self.stream_btn.setText("Stop Streaming")
            self.connect_btn.setEnabled(False)
            self.update_timer.start(10)  # Update every 10ms
            self.statusBar().showMessage('Streaming data...')
        else:
            # Stop streaming
            self.is_streaming = False
            self.stream_btn.setText("Start Streaming")
            self.connect_btn.setEnabled(True)
            self.update_timer.stop()
            self.statusBar().showMessage('Streaming stopped')

    def update_graphs(self):
        """Read data from serial and update graphs"""
        data, error = self.backend.read_line()

        if data:
            # Update all graphs with new data
            for graph in self.graphs:
                graph.update_data(data)
        elif error and "Empty line" not in error:
            self.statusBar().showMessage(f'Error: {error}')

    def clear_graphs(self):
        """Clear all graph data"""
        for graph in self.graphs:
            graph.clear_data()
        self.statusBar().showMessage('Graphs cleared')

    def closeEvent(self, event):
        """Handle window close event"""
        if self.is_streaming:
            self.toggle_streaming()

        if self.backend.serial_conn and self.backend.serial_conn.is_open:
            self.backend.disconnect()

        event.accept()


def main():
    app = QApplication(sys.argv)

    # Set dark theme for better visibility
    app.setStyle('Fusion')
    pg.setConfigOptions(antialias=True)

    dashboard = SensorDashboard(sensor_type='LPS22HHTR')
    dashboard.show()

    sys.exit(app.exec())


if __name__ == '__main__':
    main()