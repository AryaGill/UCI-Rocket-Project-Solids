from PyQt6.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt


class RPYGraph(QWidget):
    """
    Reusable widget for displaying Roll, Pitch, Yaw data with modern styling.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        self.time_data = []
        self.roll_data = []
        self.pitch_data = []
        self.yaw_data = []

        plt.style.use('dark_background')

        self.figure = Figure(figsize=(8, 5), facecolor='#1e1e1e')
        self.canvas = FigureCanvas(self.figure)
        self.axes = self.figure.add_subplot(111)
        self.figure.subplots_adjust(left=0.2, bottom=0.15)
        self.axes.set_facecolor('#2d2d2d')

        self.axes.set_title('Roll / Pitch / Yaw', fontsize=16, fontweight='bold',
                            color='#a855f7', pad=15)
        self.axes.set_xlabel('Time (s)', fontsize=12, color='#b0b0b0')
        self.axes.set_ylabel('Angle (°)', fontsize=12, color='#b0b0b0')
        self.axes.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)

        self.line_roll, = self.axes.plot([], [], color='#ff6b35', linewidth=2.5,
                                         label='Roll', antialiased=True)
        self.line_pitch, = self.axes.plot([], [], color='#28dc5e', linewidth=2.5,
                                          label='Pitch', antialiased=True)
        self.line_yaw, = self.axes.plot([], [], color='#357fff', linewidth=2.5,
                                        label='Yaw', antialiased=True)

        legend = self.axes.legend(facecolor='#2d2d2d', edgecolor='#a855f7', fontsize=10)
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
                      edgecolor='#a855f7', alpha=0.85),
            color='#e0e0e0'
        )

        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.setContentsMargins(5, 5, 5, 5)
        self.setLayout(layout)

    def _update_stats(self):
        if not self.roll_data:
            self._stats_text.set_text('')
            return
        self._stats_text.set_text(
            f"Roll  now:{self.roll_data[-1]:+7.2f}°\n"
            f"Pitch now:{self.pitch_data[-1]:+7.2f}°\n"
            f"Yaw   now:{self.yaw_data[-1]:+7.2f}°"
        )

    def update_data(self, t, roll, pitch, yaw, max_points=100):
        self.time_data.append(t)
        self.roll_data.append(roll)
        self.pitch_data.append(pitch)
        self.yaw_data.append(yaw)

        self.time_data  = self.time_data[-max_points:]
        self.roll_data  = self.roll_data[-max_points:]
        self.pitch_data = self.pitch_data[-max_points:]
        self.yaw_data   = self.yaw_data[-max_points:]

        self.line_roll.set_xdata(self.time_data);  self.line_roll.set_ydata(self.roll_data)
        self.line_pitch.set_xdata(self.time_data); self.line_pitch.set_ydata(self.pitch_data)
        self.line_yaw.set_xdata(self.time_data);   self.line_yaw.set_ydata(self.yaw_data)

        self._update_stats()
        self.axes.relim()
        self.axes.autoscale_view(True, True, True)
        self.canvas.draw()

    def clear_data(self):
        self.time_data.clear()
        self.roll_data.clear()
        self.pitch_data.clear()
        self.yaw_data.clear()

        self.line_roll.set_data([], [])
        self.line_pitch.set_data([], [])
        self.line_yaw.set_data([], [])
        self._stats_text.set_text('')

        self.axes.relim()
        self.axes.autoscale_view(True, True, True)
        self.canvas.draw_idle()