from PyQt6.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt


class AltitudeGraph(QWidget):
    """
    Reusable widget for displaying altitude data from two sensors with modern styling.
    """
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        self.time_data = []
        self.altitude_data = []
        self.filtered_altitude_data = []
        
        # Use dark background style for matplotlib
        plt.style.use('dark_background')
        
        # Create matplotlib figure with dark theme
        self.figure = Figure(figsize=(8, 5), facecolor='#1e1e1e')
        self.canvas = FigureCanvas(self.figure)
        self.axes = self.figure.add_subplot(111)
        self.figure.subplots_adjust(left=0.2, bottom=0.15)
        self.axes.set_facecolor('#2d2d2d')
        
        # Configure the plot with modern styling
        self.axes.set_title('Altitude', fontsize=16, fontweight='bold', 
                           color='#00d4ff', pad=15)
        self.axes.set_xlabel('Time (s)', fontsize=12, color='#b0b0b0')
        self.axes.set_ylabel('Altitude (m)', fontsize=12, color='#b0b0b0')
        self.axes.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
        
        # Two line plots with different colors
        self.line_raw, = self.axes.plot([], [], color='#00d4ff', linewidth=2.5, 
                                        label='Raw Alt', antialiased=True)
        self.line_filtered, = self.axes.plot([], [], color='#ff6b35', linewidth=2.5, 
                                             label='Filtered Alt', antialiased=True)
        
        # Style the legend
        legend = self.axes.legend(facecolor='#2d2d2d', edgecolor='#00d4ff', 
                                 fontsize=10)
        for text in legend.get_texts():
            text.set_color('#b0b0b0')
        
        # Style tick labels
        self.axes.tick_params(colors='#b0b0b0', labelsize=10)
        
        # Add subtle border
        for spine in self.axes.spines.values():
            spine.set_edgecolor('#404040')
            spine.set_linewidth(1)
        
        # Layout
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.setContentsMargins(5, 5, 5, 5)
        self.setLayout(layout)
        
    def update_data(self, time_value, altitude_value, filtered_altitude_value, max_points=100):
        """Update the graph with new data points from both sensors."""
        self.time_data.append(time_value)
        self.altitude_data.append(altitude_value)
        self.filtered_altitude_data.append(filtered_altitude_value)

        self.time_data = self.time_data[-max_points:]
        self.altitude_data = self.altitude_data[-max_points:]
        self.filtered_altitude_data = self.filtered_altitude_data[-max_points:]
        
        # Update raw altitude line
        self.line_raw.set_xdata(self.time_data)
        self.line_raw.set_ydata(self.altitude_data)
        
        # Update filtered altitude line
        self.line_filtered.set_xdata(self.time_data)
        self.line_filtered.set_ydata(self.filtered_altitude_data)
        
        self.axes.relim()
        self.axes.autoscale_view(True, True, True)
        
        self.canvas.draw()
        
    def clear_data(self):
        """Clear all data from the graph."""
        self.time_data = []
        self.altitude_data = []
        self.filtered_altitude_data = []
        self.line_raw.set_data([], [])
        self.line_filtered.set_data([], [])
        self.axes.relim()
        self.axes.autoscale_view(True, True, True)
        self.canvas.draw_idle()
