from PyQt6.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt


class MagGraph(QWidget):
    """
    Reusable widget for displaying temperature data with modern styling.
    """
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        self.time_data = []
        self.mag_x = []
        self.mag_y = []
        self.mag_z = []
        
        # Use dark background style
        plt.style.use('dark_background')
        
        # Create matplotlib figure with dark theme
        self.figure = Figure(figsize=(8, 5), facecolor='#1e1e1e')
        self.canvas = FigureCanvas(self.figure)
        self.axes = self.figure.add_subplot(111)
        self.figure.subplots_adjust(left=0.2, bottom=0.15)
        self.axes.set_facecolor('#2d2d2d')
        
        # Configure the plot with modern styling
        self.axes.set_title('Magnetic Field', fontsize=16, fontweight='bold', 
                           color="#38b0fb", pad=15)
        self.axes.set_xlabel('Time (s)', fontsize=12, color='#b0b0b0')
        self.axes.set_ylabel('Magnetic Field (uT)', fontsize=12, color='#b0b0b0')
        self.axes.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
        
        # Stylish line plot with warm color for temperature
        self.line_x, = self.axes.plot([], [], color='#ff6b35', linewidth=2.5, 
                                    label='X', antialiased=True)
        self.line_y, = self.axes.plot([], [], color="#28dc5e", linewidth=2.5, 
                                    label='Y', antialiased=True)
        self.line_z, = self.axes.plot([], [], color="#357fff", linewidth=2.5, 
                                    label='Z', antialiased=True)
        
        # Style the legend
        legend = self.axes.legend(facecolor='#2d2d2d', edgecolor='#38b0fb', 
                                 fontsize=10)
        legend.get_texts()[0].set_color('#b0b0b0')
        
        # Style tick labels
        self.axes.tick_params(colors='#b0b0b0', labelsize=10)
        
        # Add subtle border
        for spine in self.axes.spines.values():
            spine.set_edgecolor('#404040')
            spine.set_linewidth(1)
        
        #Current
        self._stats_text = self.axes.text(
            0.02, 0.97, '',
            transform=self.axes.transAxes,
            fontsize=9, verticalalignment='top',
            fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#1e1e1e',
                        edgecolor='#38b0fb', alpha=0.85),
            color='#e0e0e0'
        )

        # Layout
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.setContentsMargins(5, 5, 5, 5)
        self.setLayout(layout)
    
    def _update_stats(self):
        if not self.mag_x:
            self._stats_text.set_text('')
            return
        cx, cy, cz = self.mag_x[-1], self.mag_y[-1], self.mag_z[-1]
        self._stats_text.set_text(
            f"X  now:{cx:+7.2f}\n"
            f"Y  now:{cy:+7.2f}\n"
            f"Z  now:{cz:+7.2f}"
        )
        
    def update_data(self, t, x, y, z):
        """Update the graph with new data point."""
        self.time_data.append(t)
        self.mag_x.append(x)
        self.mag_y.append(y)
        self.mag_z.append(z)

        self.line_x.set_xdata(self.time_data); self.line_x.set_ydata(self.mag_x)
        self.line_y.set_xdata(self.time_data); self.line_y.set_ydata(self.mag_y)
        self.line_z.set_xdata(self.time_data); self.line_z.set_ydata(self.mag_z)

        self._update_stats()

        self.axes.relim()
        self.axes.autoscale_view(True, True, True)
        self.canvas.draw()
        
    def clear_data(self):
        """Clear all buffers and reset plot."""
        self.time_data.clear()
        self.mag_x.clear()
        self.mag_y.clear()
        self.mag_z.clear()

        self.line_x.set_data([], [])
        self.line_y.set_data([], [])
        self.line_z.set_data([], [])
        self._stats_text.set_text('')

        self.axes.relim()
        self.axes.autoscale_view(True, True, True)
        self.canvas.draw_idle()
