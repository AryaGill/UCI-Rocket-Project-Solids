import numpy as np
from PyQt6 import QtWidgets, QtCore

from motor import detect_burn_start, detect_burn_end

LAND_ALT_WINDOW = 10.0
LAND_CONFIRM    = 20


def compute_flight_events(time_ms, accel_world_z=None, altitude=None):
    events = {}

    if accel_world_z is not None:
        i_launch,  t_launch_ms  = detect_burn_start(time_ms, accel_world_z)
        i_burnout, t_burnout_ms = detect_burn_end(time_ms, accel_world_z, i_launch)
        events["launch_ms"] = float(t_launch_ms) if t_launch_ms is not None else None
        if t_burnout_ms is not None and t_launch_ms is not None:
            events["burn_time"] = float((t_burnout_ms - t_launch_ms) / 1000.0)
        else:
            events["burn_time"] = None

    if altitude is not None:
        events["max_alt"]   = float(np.max(altitude))
        apogee_idx          = int(np.argmax(altitude))
        launch_ms           = events.get("launch_ms")
        # apogee time in raw ms
        events["apogee_ms"] = float(time_ms[apogee_idx])
        events["apogee_t"]  = (events["apogee_ms"] - launch_ms) / 1000.0 if launch_ms is not None else None

        launch_alt  = float(altitude[0])
        consecutive = 0
        land_idx    = None
        for i in range(apogee_idx, len(altitude)):
            if abs(altitude[i] - launch_alt) <= LAND_ALT_WINDOW:
                consecutive += 1
                if consecutive >= LAND_CONFIRM:
                    land_idx = i - LAND_CONFIRM + 1
                    break
            else:
                consecutive = 0

        if land_idx is not None and launch_ms is not None:
            events["land_ms"]    = float(time_ms[land_idx])
            events["flight_time"] = (events["land_ms"] - launch_ms) / 1000.0
        else:
            events["land_ms"]    = None
            events["flight_time"] = None

    return events


class FlightEventsBar(QtWidgets.QWidget):
    markers_toggled = QtCore.pyqtSignal(bool, dict)

    _TILE_STYLE = (
        "background:#080C14; border:1px solid #263245; border-radius:7px; "
        "padding:8px 14px;"
    )
    _VAL_STYLE  = (
        "font-family:'Menlo','Monaco',monospace; font-size:16px; "
        "font-weight:800; color:{color};"
    )
    _LBL_STYLE  = "font-size:10px; font-weight:650; color:#91A0B5;"

    _FIELDS = [
        ("burn_time",   "Burn Time",    "#F59E0B"),
        ("max_alt",     "Max Altitude", "#22C55E"),
        ("flight_time", "Flight Time",  "#38BDF8"),
    ]

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("Panel")
        self._events = {}
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(12, 12, 12, 12)
        outer.setSpacing(10)

        header = QtWidgets.QHBoxLayout()
        title = QtWidgets.QLabel("Flight Events")
        title.setObjectName("SectionTitle")
        header.addWidget(title)
        header.addStretch()
        self.display_btn = QtWidgets.QPushButton("Graph Markers")
        self.display_btn.setCheckable(True)
        self.display_btn.setFixedHeight(30)
        self.display_btn.setFixedWidth(130)
        self.display_btn.toggled.connect(self._on_toggle)
        header.addWidget(self.display_btn)
        outer.addLayout(header)

        tile_row = QtWidgets.QHBoxLayout()
        tile_row.setSpacing(8)
        self._tiles = {}
        for key, label, color in self._FIELDS:
            tile = QtWidgets.QWidget()
            tile.setStyleSheet(self._TILE_STYLE)
            vbox = QtWidgets.QVBoxLayout(tile)
            vbox.setContentsMargins(0, 0, 0, 0)
            vbox.setSpacing(1)
            val_lbl = QtWidgets.QLabel("--")
            val_lbl.setStyleSheet(self._VAL_STYLE.format(color=color))
            val_lbl.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
            cap_lbl = QtWidgets.QLabel(label)
            cap_lbl.setStyleSheet(self._LBL_STYLE)
            cap_lbl.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
            vbox.addWidget(val_lbl)
            vbox.addWidget(cap_lbl)
            tile_row.addWidget(tile, stretch=1)
            self._tiles[key] = val_lbl
        outer.addLayout(tile_row)
        self.setFixedHeight(128)

    def _on_toggle(self, checked):
        self.markers_toggled.emit(checked, self._events)

    def update_events(self, events):
        self._events = events
        def fmt_t(v): return ("T+%.2f s" % v) if v is not None else "--"
        def fmt_m(v): return ("%.1f m"   % v) if v is not None else "--"
        self._tiles["burn_time"].setText(fmt_t(events.get("burn_time")))
        self._tiles["max_alt"].setText(fmt_m(events.get("max_alt")))
        self._tiles["flight_time"].setText(fmt_t(events.get("flight_time")))
        if self.display_btn.isChecked():
            self.markers_toggled.emit(True, self._events)

    def clear(self):
        self._events = {}
        for lbl in self._tiles.values():
            lbl.setText("--")
        if self.display_btn.isChecked():
            self.display_btn.setChecked(False)
