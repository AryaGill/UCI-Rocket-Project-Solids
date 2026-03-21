from PyQt6.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt


class APDeploymentGraph(QWidget):
    """
    Reusable widget for displaying AP deployment percentage
    and predicted apogee vs time with modern styling.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        self.time_data = []
        self.ap_data = []
        self.pred_apo_data = []

        plt.style.use('dark_background')

        self.figure = Figure(figsize=(8, 5), facecolor='#1e1e1e')
        self.canvas = FigureCanvas(self.figure)

        # Two y-axes: left for AP deployment, right for predicted apogee
        self.axes = self.figure.add_subplot(111)
        self.axes_apo = self.axes.twinx()

        self.figure.subplots_adjust(left=0.15, bottom=0.15, right=0.85)
        self.axes.set_facecolor('#2d2d2d')

        # ── Left axis: AP Deployment ──────────────────────────────────────────
        self.axes.set_title('AP Deployment & Predicted Apogee vs Time',
                            fontsize=14, fontweight='bold',
                            color='#a855f7', pad=15)
        self.axes.set_xlabel('Time (s)', fontsize=12, color='#b0b0b0')
        self.axes.set_ylabel('AP Deployment (%)', fontsize=12, color='#a855f7')
        self.axes.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
        self.axes.tick_params(axis='y', colors='#a855f7', labelsize=10)
        self.axes.tick_params(axis='x', colors='#b0b0b0', labelsize=10)

        self.line_ap, = self.axes.plot([], [], color='#a855f7', linewidth=2.5,
                                       label='AP Deployment', antialiased=True)

        # ── Right axis: Predicted Apogee ─────────────────────────────────────
        self.axes_apo.set_ylabel('Predicted Apogee (m)', fontsize=12, color='#00d4ff')
        self.axes_apo.tick_params(axis='y', colors='#00d4ff', labelsize=10)

        self.line_apo, = self.axes_apo.plot([], [], color='#00d4ff', linewidth=2.0,
                                             linestyle='--', label='Pred Apogee',
                                             antialiased=True)

        # ── Spines ───────────────────────────────────────────────────────────
        for spine in self.axes.spines.values():
            spine.set_edgecolor('#404040')
            spine.set_linewidth(1)
        for spine in self.axes_apo.spines.values():
            spine.set_edgecolor('#404040')
            spine.set_linewidth(1)

        # ── Combined legend ──────────────────────────────────────────────────
        lines = [self.line_ap, self.line_apo]
        labels = [l.get_label() for l in lines]
        legend = self.axes.legend(lines, labels,
                                  facecolor='#2d2d2d', edgecolor='#a855f7',
                                  fontsize=10, loc='upper left')
        for text in legend.get_texts():
            text.set_color('#b0b0b0')

        # ── Stats annotation ─────────────────────────────────────────────────
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

    # ── Internals ─────────────────────────────────────────────────────────────

    def _update_stats(self):
        if not self.ap_data:
            self._stats_text.set_text('')
            return
        ap_now = self.ap_data[-1]
        apo_now = self.pred_apo_data[-1] if self.pred_apo_data else None

        lines = [f"AP Deploy: {ap_now:+.2f} %"]
        if apo_now is not None:
            lines.append(f"Pred Apo:  {apo_now:+.1f} m")
        self._stats_text.set_text('\n'.join(lines))

    # ── Public API ────────────────────────────────────────────────────────────

    def update_data(self, t, ap_deployment, pred_apo=None, max_points=100):
        """
        Update the graph with a new data point.

        Args:
            t:             Time value (seconds).
            ap_deployment: AP deployment value.
            pred_apo:      Predicted apogee in metres (optional).
            max_points:    Rolling window size.
        """
        self.time_data.append(t)
        self.ap_data.append(ap_deployment)
        self.pred_apo_data.append(pred_apo)

        self.time_data     = self.time_data[-max_points:]
        self.ap_data       = self.ap_data[-max_points:]
        self.pred_apo_data = self.pred_apo_data[-max_points:]

        self.line_ap.set_xdata(self.time_data)
        self.line_ap.set_ydata(self.ap_data)

        # Only plot apogee points where data is not None
        valid_apo = [(t, a) for t, a in zip(self.time_data, self.pred_apo_data)
                     if a is not None]
        if valid_apo:
            t_vals, a_vals = zip(*valid_apo)
            self.line_apo.set_xdata(list(t_vals))
            self.line_apo.set_ydata(list(a_vals))
        else:
            self.line_apo.set_xdata([])
            self.line_apo.set_ydata([])

        self._update_stats()

        self.axes.relim()
        self.axes.autoscale_view(True, True, True)
        self.axes_apo.relim()
        self.axes_apo.autoscale_view(True, True, True)
        self.canvas.draw()

    def clear_data(self):
        """Clear all buffers and reset the plot."""
        self.time_data.clear()
        self.ap_data.clear()
        self.pred_apo_data.clear()

        self.line_ap.set_data([], [])
        self.line_apo.set_data([], [])
        self._stats_text.set_text('')

        self.axes.relim()
        self.axes.autoscale_view(True, True, True)
        self.axes_apo.relim()
        self.axes_apo.autoscale_view(True, True, True)
        self.canvas.draw_idle()