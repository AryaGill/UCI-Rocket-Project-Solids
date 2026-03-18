from PyQt6.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, QLabel, QFrame
from PyQt6.QtCore import Qt


# ADC parameters
ADC_MAX_COUNTS = 65536       # 16-bit ADC full scale
ADC_REF_VOLTAGE = 3.3        # volts at full scale

# Raw ADC count above this is considered "connected" (ematch present)
# Corresponds to ~2.0 V: (2.0 / 3.3) * 65536 ≈ 39,759
CONNECTED_THRESHOLD_COUNTS = int((2.0 / ADC_REF_VOLTAGE) * ADC_MAX_COUNTS)


def adc_to_volts(raw: float) -> float:
    """Convert a raw 16-bit ADC count to volts."""
    return (raw / ADC_MAX_COUNTS) * ADC_REF_VOLTAGE


class _EMatchIndicator(QWidget):
    """
    Single ematch channel indicator.
    Shows channel label, live voltage reading, and a colour-coded
    connectivity dot (green = connected, red = open / no ematch).
    """

    def __init__(self, label: str, color: str, parent=None):
        super().__init__(parent)
        self._accent = color

        root = QVBoxLayout(self)
        root.setContentsMargins(8, 6, 8, 6)
        root.setSpacing(4)

        # ── channel name ────────────────────────────────────────────────────
        name_lbl = QLabel(label)
        name_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        name_lbl.setStyleSheet(f"color: {color}; font-weight: bold; font-size: 11px;")
        root.addWidget(name_lbl)

        # ── voltage reading ──────────────────────────────────────────────────
        self.voltage_lbl = QLabel("-.-- V")
        self.voltage_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.voltage_lbl.setStyleSheet(
            "color: #e0e0e0; font-size: 13px; font-weight: bold; font-family: monospace;"
        )
        root.addWidget(self.voltage_lbl)

        # ── status dot + text ────────────────────────────────────────────────
        status_row = QHBoxLayout()
        status_row.setSpacing(5)
        status_row.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.dot = QFrame()
        self.dot.setFixedSize(12, 12)
        self.dot.setFrameShape(QFrame.Shape.NoFrame)
        status_row.addWidget(self.dot)

        self.status_lbl = QLabel("No Data")
        self.status_lbl.setStyleSheet("font-size: 10px; color: #808080;")
        status_row.addWidget(self.status_lbl)

        root.addLayout(status_row)

        # start in the "unknown" state
        self._set_unknown()

        # outer card border
        self.setStyleSheet(f"""
            _EMatchIndicator, QWidget {{
                background-color: #2a2a2a;
            }}
        """)
        self.setAutoFillBackground(True)

    # ── internal state helpers ───────────────────────────────────────────────

    def _set_dot_color(self, hex_color: str):
        self.dot.setStyleSheet(
            f"border-radius: 6px; background-color: {hex_color};"
        )

    def _set_unknown(self):
        self._set_dot_color("#606060")
        self.status_lbl.setText("No Data")
        self.status_lbl.setStyleSheet("font-size: 10px; color: #808080;")

    # ── public API ───────────────────────────────────────────────────────────

    def update_voltage(self, raw: float):
        """Refresh the indicator with a raw ADC count (0–65536)."""
        voltage = adc_to_volts(raw)
        self.voltage_lbl.setText(f"{voltage:.2f} V")

        if raw >= CONNECTED_THRESHOLD_COUNTS:
            self._set_dot_color("#28dc5e")          # green  – connected
            self.status_lbl.setText("Connected")
            self.status_lbl.setStyleSheet("font-size: 10px; color: #28dc5e;")
        else:
            self._set_dot_color("#ff4040")          # red    – open circuit
            self.status_lbl.setText("Open / No Ematch")
            self.status_lbl.setStyleSheet("font-size: 10px; color: #ff4040;")


class EMatchPanel(QWidget):
    """
    Horizontal bar showing one indicator per pyro channel.

    Channels (matching backend COLUMNS names):
        main_p_ematch_voltage
        main_s_ematch_voltage
        drogue_p_ematch_voltage
        drogue_s_ematch_voltage
    """

    # Maps data-dict key  →  (display label,  accent colour)
    CHANNELS = [
        ("drogue_p_ematch_voltage", "Drogue\nPrimary",   "#ff6b35"),
        ("drogue_s_ematch_voltage", "Drogue\nSecondary", "#ff9966"),
        ("main_p_ematch_voltage",   "Main\nPrimary",     "#28dc5e"),
        ("main_s_ematch_voltage",   "Main\nSecondary",   "#66ee99"),
    ]

    def __init__(self, parent=None):
        super().__init__(parent)

        outer = QHBoxLayout(self)
        outer.setContentsMargins(6, 4, 6, 4)
        outer.setSpacing(0)

        # # ── section title ────────────────────────────────────────────────────
        # title = QLabel("E-Match\nStatus")
        # title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        # title.setStyleSheet(
        #     "color: #b0b0b0; font-size: 10px; font-weight: bold; padding: 0 8px;"
        # )
        # outer.addWidget(title)

        # thin divider
        div = QFrame()
        div.setFrameShape(QFrame.Shape.VLine)
        div.setStyleSheet("color: #404040;")
        outer.addWidget(div)

        # ── one card per channel ─────────────────────────────────────────────
        self._indicators: dict[str, _EMatchIndicator] = {}

        for key, label, color in self.CHANNELS:
            ind = _EMatchIndicator(label, color)
            outer.addWidget(ind)
            self._indicators[key] = ind

            if key != self.CHANNELS[-1][0]:
                sep = QFrame()
                sep.setFrameShape(QFrame.Shape.VLine)
                sep.setFixedWidth(1)
                sep.setStyleSheet("background-color: #404040;")
                outer.addWidget(sep)

        self.setStyleSheet("background-color: #232323; border-radius: 6px;")
        self.setFixedHeight(90)
        self.setMaximumWidth(400)

    # ── public API ───────────────────────────────────────────────────────────

    def update_data(self, data: dict):
        """
        Feed a telemetry data dict.  Only keys that are present and
        numeric will update the corresponding indicator.
        """
        for key, _label, _color in self.CHANNELS:
            value = data.get(key)
            if value is not None:
                try:
                    self._indicators[key].update_voltage(float(value))
                except (ValueError, TypeError):
                    pass