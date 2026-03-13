from PyQt6.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt


class VelocityGraph(QWidget):
    """
    Reusable widget for displaying world-frame velocity data with modern styling.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        self.time_data = []
        self.vel_x = []
        self.vel_y = []
        self.vel_z = []

        plt.style.use('dark_background')

        self.figure = Figure(figsize=(8, 5), facecolor='#1e1e1e')
        self.canvas = FigureCanvas(self.figure)
        self.axes = self.figure.add_subplot(111)
        self.figure.subplots_adjust(left=0.2, bottom=0.15)
        self.axes.set_facecolor('#2d2d2d')

        self.axes.set_title('Velocity', fontsize=16, fontweight='bold',
                            color='#00d4ff', pad=15)
        self.axes.set_xlabel('Time (s)', fontsize=12, color='#b0b0b0')
        self.axes.set_ylabel('Velocity (m/s)', fontsize=12, color='#b0b0b0')
        self.axes.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)

        self.line_x, = self.axes.plot([], [], color='#ff6b35', linewidth=2.5,
                                      label='X', antialiased=True)
        self.line_y, = self.axes.plot([], [], color='#28dc5e', linewidth=2.5,
                                      label='Y', antialiased=True)
        self.line_z, = self.axes.plot([], [], color='#357fff', linewidth=2.5,
                                      label='Z', antialiased=True)

        legend = self.axes.legend(facecolor='#2d2d2d', edgecolor='#00d4ff', fontsize=10)
        for text in legend.get_texts():
            text.set_color('#b0b0b0')

        self.axes.tick_params(colors='#b0b0b0', labelsize=10)
        for spine in self.axes.spines.values():
            spine.set_edgecolor('#404040')
            spine.set_linewidth(1)

        self._stats_text = self.axes.text(
            0.02, 0.97, '',
            transform=self.axes.transAxes,
            fontsize=9, verticalalignment='top',
            fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#1e1e1e',
                      edgecolor='#00d4ff', alpha=0.85),
            color='#e0e0e0'
        )

        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.setContentsMargins(5, 5, 5, 5)
        self.setLayout(layout)

    def _update_stats(self):
        if not self.vel_x:
            self._stats_text.set_text('')
            return
        self._stats_text.set_text(
            f"X now:{self.vel_x[-1]:+7.2f} m/s\n"
            f"Y now:{self.vel_y[-1]:+7.2f} m/s\n"
            f"Z now:{self.vel_z[-1]:+7.2f} m/s"
        )

    def update_data(self, t, vx, vy, vz, max_points=100):
        self.time_data.append(t)
        self.vel_x.append(vx)
        self.vel_y.append(vy)
        self.vel_z.append(vz)

        self.time_data = self.time_data[-max_points:]
        self.vel_x     = self.vel_x[-max_points:]
        self.vel_y     = self.vel_y[-max_points:]
        self.vel_z     = self.vel_z[-max_points:]

        self.line_x.set_xdata(self.time_data); self.line_x.set_ydata(self.vel_x)
        self.line_y.set_xdata(self.time_data); self.line_y.set_ydata(self.vel_y)
        self.line_z.set_xdata(self.time_data); self.line_z.set_ydata(self.vel_z)

        self._update_stats()
        self.axes.relim()
        self.axes.autoscale_view(True, True, True)
        self.canvas.draw()

    def clear_data(self):
        self.time_data.clear()
        self.vel_x.clear()
        self.vel_y.clear()
        self.vel_z.clear()

        self.line_x.set_data([], [])
        self.line_y.set_data([], [])
        self.line_z.set_data([], [])
        self._stats_text.set_text('')

        self.axes.relim()
        self.axes.autoscale_view(True, True, True)
        self.canvas.draw_idle()