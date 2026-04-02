from PyQt6.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt


class APGraph(QWidget):
    """
    Reusable widget for displaying apogee data with modern styling.
    """
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        self.time_data = []
        self.ap_data = []
        # Use dark background style
        plt.style.use('dark_background')
        
        # Create matplotlib figure with dark theme
        self.figure = Figure(figsize=(8, 5), facecolor='#1e1e1e')
        self.canvas = FigureCanvas(self.figure)
        self.axes = self.figure.add_subplot(111)
        self.figure.subplots_adjust(left=0.2, bottom=0.15)
        self.axes.set_facecolor('#2d2d2d')
        
        # Configure the plot with modern styling
        self.axes.set_title('Predicted Apogee', fontsize=16, fontweight='bold', 
                           color="#7afdff", pad=15)
        self.axes.set_xlabel('Time (s)', fontsize=12, color='#b0b0b0')
        self.axes.set_ylabel('Predicted Apogee (ft)', fontsize=12, color='#b0b0b0')
        self.axes.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
        
        # Stylish line plot with warm color for Apogee
        self.line, = self.axes.plot([], [], color='#7afdff', linewidth=2.5, 
                                    label='Predicted Apogee', antialiased=True)
        
        # Style the legend
        legend = self.axes.legend(facecolor='#2d2d2d', edgecolor='#7afdff', 
                                 fontsize=10)
        legend.get_texts()[0].set_color('#b0b0b0')
        
        # Style tick labels
        self.axes.tick_params(colors='#b0b0b0', labelsize=10)
        
        # Add subtle border
        for spine in self.axes.spines.values():
            spine.set_edgecolor('#404040')
            spine.set_linewidth(1)

        # Current / max value annotation (top-left, inside axes)
        self._stats_text = self.axes.text(
            0.02, 0.97, '',
            transform=self.axes.transAxes,
            fontsize=9, verticalalignment='top',
            fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#1e1e1e',
                      edgecolor='#7afdff', alpha=0.85),
            color='#e0e0e0'
        )
        
        # Layout
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.setContentsMargins(5, 5, 5, 5)
        self.setLayout(layout)

    def _update_stats(self):
        if not self.ap_data:
            self._stats_text.set_text('')
            return
        self._stats_text.set_text(f"Now: {self.ap_data[-1]:+.2f} ft")
        
    def update_data(self, time_value, ap_value, max_points=100):
        """Update the graph with new data point."""
        self.time_data.append(time_value)
        self.ap_data.append(ap_value)

        self.time_data = self.time_data[-max_points:]
        self.ap_data = self.ap_data[-max_points:]
        
        self.line.set_xdata(self.time_data)
        self.line.set_ydata(self.ap_data)
        
        self._update_stats()

        self.axes.relim()
        self.axes.autoscale_view(True, True, True)
        
        self.canvas.draw()
        
    def clear_data(self):
        """Clear all data from the graph."""
        self.time_data = []
        self.ap_data = []
        self.line.set_xdata([])
        self.line.set_ydata([])
        self._stats_text.set_text('')
        self.canvas.draw()