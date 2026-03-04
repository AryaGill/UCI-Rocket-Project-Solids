from PyQt6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QLabel, QGridLayout, QFrame
from PyQt6.QtCore import Qt


class MaxValuesTable(QWidget):
    """
    Compact panel showing the all-time maximum values recorded across
    every active graph during the current session.

    Layout (one row per sensor group, XYZ columns where applicable):

        Sensor          |  Value / X       |  Y         |  Z
        ────────────────┼──────────────────┼────────────┼──────────
        Alt Raw         |  0.00 m          |            |
        Alt Filtered    |  0.00 m          |            |
        Temperature     |  0.00 °C         |            |
        Accel LIS       |  0.00 m/s²       | 0.00 m/s²  | 0.00 m/s²
        Accel World     |  0.00 m/s²       | 0.00 m/s²  | 0.00 m/s²
        Angular Speed   |  0.00 rad/s      | 0.00 rad/s | 0.00 rad/s
    """

    # ── colour accents matching each graph ──────────────────────────────────
    _ACCENT = {
        "alt":         "#00d4ff",
        "temp":        "#ff6b35",
        "accel_lis":   "#00b71f",
        "accel_world": "#00b71f",
        "angular":     "#fffb00",
    }

    def __init__(self, parent=None):
        super().__init__(parent)
        self._build_ui()
        self.reset()

    # ── UI construction ──────────────────────────────────────────────────────

    def _build_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 2, 4, 2)
        outer.setSpacing(4)

        # Title
        title = QLabel("Peak Values")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title.setFixedHeight(30)
        title.setStyleSheet(
            "color: #b0b0b0; font-size: 10px; font-weight: bold;"
        )
        outer.addWidget(title)

        # Grid table
        grid = QGridLayout()
        grid.setHorizontalSpacing(12)
        grid.setVerticalSpacing(2)
        outer.addLayout(grid)

        # Column headers
        for col, text in enumerate(["Sensor", "Value / X", "Y", "Z"]):
            lbl = QLabel(text)
            lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            lbl.setStyleSheet(
                "color: #606060; font-size: 9px; font-weight: bold;"
            )
            grid.addWidget(lbl, 0, col)

        # Thin separator under header
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet("color: #404040;")
        grid.addWidget(sep, 1, 0, 1, 4)

        # Row definitions: (attr_key, display_label, accent, unit, has_xyz)
        rows = [
            ("alt_raw",      "Alt Raw",      "alt",         "m",      False),
            ("alt_filt",     "Alt Filtered", "alt",         "m",      False),
            ("temp",         "Temperature",  "temp",        "°C",     False),
            ("accel_lis",    "Accel LIS",    "accel_lis",   "m/s²",   True),
            ("accel_world",  "Accel World",  "accel_world", "m/s²",   True),
            ("angular",      "Angular Spd",  "angular",     "rad/s",  True),
        ]

        self._cells: dict[str, list[QLabel]] = {}   # key → [val_lbl] or [x_lbl, y_lbl, z_lbl]

        for i, (key, label, accent, unit, has_xyz) in enumerate(rows):
            row = i + 2          # rows 0=header, 1=sep, 2+= data
            color = self._ACCENT[accent]

            # Sensor name
            name_lbl = QLabel(label)
            name_lbl.setStyleSheet(
                f"color: {color}; font-size: 9px; font-weight: bold;"
            )
            grid.addWidget(name_lbl, row, 0)

            if has_xyz:
                cell_lbls = []
                for col in range(1, 4):
                    lbl = QLabel(f"-.-- {unit}")
                    lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
                    lbl.setStyleSheet(
                        "color: #c0c0c0; font-size: 9px; font-family: monospace;"
                    )
                    lbl.setMinimumWidth(72)
                    grid.addWidget(lbl, row, col)
                    cell_lbls.append(lbl)
                self._cells[key] = cell_lbls
            else:
                val_lbl = QLabel(f"-.-- {unit}")
                val_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
                val_lbl.setStyleSheet(
                    "color: #c0c0c0; font-size: 9px; font-family: monospace;"
                )
                val_lbl.setMinimumWidth(72)
                grid.addWidget(val_lbl, row, 1, 1, 3)   # span remaining columns
                self._cells[key] = [val_lbl]

        self.setStyleSheet("background-color: #232323; border-radius: 6px;")

    # ── public API ───────────────────────────────────────────────────────────

    def reset(self):
        """Clear all stored maxima back to None."""
        self._max: dict[str, list] = {
            "alt_raw":      [None],
            "alt_filt":     [None],
            "temp":         [None],
            "accel_lis":    [None, None, None],
            "accel_world":  [None, None, None],
            "angular":      [None, None, None],
        }
        self._units = {
            "alt_raw":      "m",
            "alt_filt":     "m",
            "temp":         "°C",
            "accel_lis":    "m/s²",
            "accel_world":  "m/s²",
            "angular":      "rad/s",
        }
        self._refresh_all()

    def update_data(self, data: dict):
        """
        Feed a telemetry data dict (same format as handle_new_data).
        Tracks maximums and refreshes only changed cells.
        """
        changed = set()

        def _maybe_update(key: str, idx: int, value):
            if value is None:
                return
            try:
                v = float(value)
            except (TypeError, ValueError):
                return
            if self._max[key][idx] is None or v > self._max[key][idx]:
                self._max[key][idx] = v
                changed.add(key)

        _maybe_update("alt_raw",   0, data.get("Alt"))
        _maybe_update("alt_filt",  0, data.get("Filtered_Alt"))
        _maybe_update("temp",      0, data.get("Temp"))

        _maybe_update("accel_lis", 0, data.get("Accel_X1"))
        _maybe_update("accel_lis", 1, data.get("Accel_Y1"))
        _maybe_update("accel_lis", 2, data.get("Accel_Z1"))

        _maybe_update("accel_world", 0, data.get("Accel_world_x"))
        _maybe_update("accel_world", 1, data.get("Accel_world_y"))
        _maybe_update("accel_world", 2, data.get("Accel_world_z"))

        _maybe_update("angular", 0, data.get("Gyro_X"))
        _maybe_update("angular", 1, data.get("Gyro_Y"))
        _maybe_update("angular", 2, data.get("Gyro_Z"))

        for key in changed:
            self._refresh_row(key)

    # ── internals ────────────────────────────────────────────────────────────

    def _fmt(self, value, unit: str) -> str:
        if value is None:
            return f"-.-- {unit}"
        return f"{value:+.2f} {unit}"

    def _refresh_row(self, key: str):
        unit = self._units[key]
        cells = self._cells[key]
        maxes = self._max[key]

        if len(cells) == 1:
            cells[0].setText(self._fmt(maxes[0], unit))
        else:
            for lbl, val in zip(cells, maxes):
                lbl.setText(self._fmt(val, unit))

    def _refresh_all(self):
        for key in self._cells:
            self._refresh_row(key)