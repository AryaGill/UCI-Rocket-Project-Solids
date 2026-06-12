import sys
import csv
import numpy as np
from PyQt6 import QtWidgets, QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl

from Frontend.motor import derive_thrust, compute_motor_stats
from Frontend.events_bar import FlightEventsBar, compute_flight_events

QtWidgets.QApplication.setAttribute(QtCore.Qt.ApplicationAttribute.AA_ShareOpenGLContexts)

fmt = QtGui.QSurfaceFormat()
fmt.setVersion(4, 1)
fmt.setProfile(QtGui.QSurfaceFormat.OpenGLContextProfile.CoreProfile)
fmt.setDepthBufferSize(24)
QtGui.QSurfaceFormat.setDefaultFormat(fmt)

COLORS = [
    "#e74c3c", "#3498db", "#2ecc71", "#f39c12", "#9b59b6",
    "#1abc9c", "#e67e22", "#34495e", "#e91e63", "#00bcd4",
    "#8bc34a", "#ff5722", "#607d8b", "#795548", "#ffeb3b",
]

def color_for_index(i):
    return COLORS[i % len(COLORS)]


def build_rocket_mesh():
    verts, faces, colors_list = [], [], []

    def add_cylinder(r, z_bot, z_top, segs, color):
        base = len(verts)
        for i in range(segs):
            a = 2 * np.pi * i / segs
            verts.append([r * np.cos(a), r * np.sin(a), z_bot])
        for i in range(segs):
            a = 2 * np.pi * i / segs
            verts.append([r * np.cos(a), r * np.sin(a), z_top])
        for i in range(segs):
            n = (i + 1) % segs
            faces.append([base + i, base + n, base + segs + n])
            faces.append([base + i, base + segs + n, base + segs + i])
            colors_list.extend([color, color])

    def add_cone(r, z_base, z_tip, segs, color):
        base = len(verts)
        for i in range(segs):
            a = 2 * np.pi * i / segs
            verts.append([r * np.cos(a), r * np.sin(a), z_base])
        tip_idx = len(verts)
        verts.append([0, 0, z_tip])
        for i in range(segs):
            n = (i + 1) % segs
            faces.append([base + i, base + n, tip_idx])
            colors_list.append(color)

    def add_fin(pts, color):
        base = len(verts)
        verts.extend(pts)
        faces.append([base, base + 1, base + 2])
        faces.append([base, base + 2, base + 3])
        colors_list.extend([color, color])

    body_r = 0.15
    body_bot, body_top, nose_tip = -1.0, 0.8, 1.5
    silver = (0.75, 0.75, 0.80, 1.0)
    red    = (0.9,  0.2,  0.2,  1.0)
    dark   = (0.3,  0.3,  0.35, 1.0)

    add_cylinder(body_r, body_bot, body_top, 24, silver)
    add_cone(body_r, body_top, nose_tip, 24, red)
    add_cone(body_r, body_bot, body_bot - 0.12, 24, dark)

    for ang_deg in [0, 90, 180, 270]:
        ang = np.radians(ang_deg)
        cx, cy = np.cos(ang), np.sin(ang)
        ox, oy = -np.sin(ang) * 0.01, np.cos(ang) * 0.01
        fw, fh = 0.35, 0.55
        inner, outer = body_r, body_r + fw
        add_fin([
            [cx * inner + ox, cy * inner + oy, body_bot],
            [cx * outer + ox, cy * outer + oy, body_bot],
            [cx * outer + ox, cy * outer + oy, body_bot + fh],
            [cx * inner + ox, cy * inner + oy, body_bot + fh],
        ], red)

    return (np.array(verts, dtype=np.float32),
            np.array(faces, dtype=np.uint32),
            np.array(colors_list, dtype=np.float32))


class MappingDialog(QtWidgets.QDialog):
    GUESSES = {
        "time":          ["time", "t", "timestamp", "time_s", "time_ms"],
        "qw":            ["quat_w", "q_w", "qw", "w", "q0"],
        "qx":            ["quat_x", "q_x", "qx", "q1"],
        "qy":            ["quat_y", "q_y", "qy", "q2"],
        "qz":            ["quat_z", "q_z", "qz", "q3"],
        "accel_world_z": ["accel_world_z", "az_world", "a_world_z", "accel_z_world"],
        "altitude":      ["altitude", "alt", "altitude_m", "baro_alt", "height"],
    }
    LABELS = {
        "time":          "Time axis",
        "qw":            "Quaternion W",
        "qx":            "Quaternion X",
        "qy":            "Quaternion Y",
        "qz":            "Quaternion Z",
        "accel_world_z": "Accel World Z (m/s^2)",
        "altitude":      "Altitude (m)",
    }
    GROUPS = [
        ("General",    ["time"]),
        ("Quaternion", ["qw", "qx", "qy", "qz"]),
        ("Motor",      ["accel_world_z"]),
        ("3D Flight",  ["altitude"]),
    ]

    def __init__(self, headers, current_mapping, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Column Mapping")
        self.setMinimumWidth(440)
        self.headers   = headers
        self.combos    = {}
        self.mass_spin = None
        layout = QtWidgets.QVBoxLayout(self)
        layout.setSpacing(12)
        hint = QtWidgets.QLabel("Map CSV columns to program fields. Leave '-- none --' to skip.")
        hint.setStyleSheet("color:#a6adc8; font-size:12px;")
        layout.addWidget(hint)
        for group_name, fields in self.GROUPS:
            box  = QtWidgets.QGroupBox(group_name)
            form = QtWidgets.QFormLayout(box)
            for field in fields:
                cb = QtWidgets.QComboBox()
                cb.addItem("-- none --")
                for h in headers:
                    cb.addItem(h)
                cur = current_mapping.get(field)
                if cur and cur in headers:
                    cb.setCurrentText(cur)
                else:
                    for guess in self.GUESSES.get(field, []):
                        match = next((h for h in headers if h.lower() == guess), None)
                        if match:
                            cb.setCurrentText(match)
                            break
                self.combos[field] = cb
                form.addRow(self.LABELS[field] + ":", cb)
            if group_name == "Motor":
                self.mass_spin = QtWidgets.QDoubleSpinBox()
                self.mass_spin.setRange(0.01, 9999.0)
                self.mass_spin.setDecimals(3)
                self.mass_spin.setSuffix("  kg")
                self.mass_spin.setValue(current_mapping.get("mass_kg", 1.0) or 1.0)
                form.addRow("Rocket Mass:", self.mass_spin)
            layout.addWidget(box)
        btns = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok |
            QtWidgets.QDialogButtonBox.StandardButton.Cancel
        )
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)

    def get_mapping(self):
        result = {}
        for field, cb in self.combos.items():
            val = cb.currentText()
            result[field] = val if val != "-- none --" else None
        if self.mass_spin:
            result["mass_kg"] = self.mass_spin.value()
        return result


class RocketView(gl.GLViewWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setCameraPosition(distance=5, elevation=20, azimuth=45)
        grid = gl.GLGridItem()
        grid.setSize(6, 6)
        grid.setSpacing(1, 1)
        grid.setColor((60, 60, 80, 120))
        self.addItem(grid)
        self.addItem(gl.GLAxisItem(size=QtGui.QVector3D(1.0, 1.0, 2.5)))
        verts, faces, face_colors = build_rocket_mesh()
        md = gl.MeshData(vertexes=verts, faces=faces, faceColors=face_colors)
        self.rocket = gl.GLMeshItem(meshdata=md, smooth=True, shader="normalColor", glOptions="opaque")
        self.addItem(self.rocket)
        self.info_label = QtWidgets.QLabel("qw=1.000  qx=0.000  qy=0.000  qz=0.000", self)
        self.info_label.setStyleSheet(
            "color:#cdd6f4; background:rgba(24,24,37,180); "
            "font-family:monospace; font-size:11px; padding:4px 8px; border-radius:4px;"
        )
        self.info_label.move(8, 8)
        self.info_label.adjustSize()

    def set_quaternion(self, w, x, y, z):
        q = QtGui.QQuaternion(w, x, y, z).normalized()
        axis, angle = q.getAxisAndAngle()
        self.rocket.resetTransform()
        self.rocket.rotate(angle, axis.x(), axis.y(), axis.z())
        self.info_label.setText("qw=%+.3f  qx=%+.3f  qy=%+.3f  qz=%+.3f" % (w, x, y, z))
        self.info_label.adjustSize()


class MotorPanel(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)
        title = QtWidgets.QLabel("Motor Classification")
        title.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        title.setStyleSheet("font-weight:bold; font-size:13px; color:#cdd6f4;")
        layout.addWidget(title)
        self.designation_label = QtWidgets.QLabel("--")
        self.designation_label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.designation_label.setStyleSheet(
            "font-size:42px; font-weight:bold; color:#f39c12; "
            "background:#181825; border-radius:8px; padding:8px;"
        )
        layout.addWidget(self.designation_label)
        grid = QtWidgets.QGridLayout()
        grid.setSpacing(6)
        self.stat_labels = {}
        for row, (key, label, unit) in enumerate([
            ("total_impulse", "Total Impulse", "N*s"),
            ("burn_time",     "Burn Time",     "s"),
            ("avg_thrust",    "Avg Thrust",    "N"),
            ("peak_thrust",   "Peak Thrust",   "N"),
        ]):
            nl = QtWidgets.QLabel(label)
            nl.setStyleSheet("color:#a6adc8; font-size:11px;")
            vl = QtWidgets.QLabel("--")
            vl.setStyleSheet(
                "color:#cdd6f4; font-family:monospace; font-size:12px; font-weight:bold;"
                "background:#181825; border-radius:4px; padding:2px 6px;"
            )
            ul = QtWidgets.QLabel(unit)
            ul.setStyleSheet("color:#585b70; font-size:11px;")
            grid.addWidget(nl, row, 0)
            grid.addWidget(vl, row, 1)
            grid.addWidget(ul, row, 2)
            self.stat_labels[key] = vl
        layout.addLayout(grid)
        layout.addStretch()

    def update_stats(self, stats):
        d = stats.get("designation")
        self.designation_label.setText(d if d else "--")
        for key, lbl in self.stat_labels.items():
            val = stats.get(key)
            lbl.setText("%.3f" % val if val is not None else "--")

    def clear(self):
        self.designation_label.setText("--")
        for lbl in self.stat_labels.values():
            lbl.setText("--")


class FlightReviewer(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Rocket Flight Reviewer")
        self.resize(1700, 860)
        self.data           = {}
        self.headers        = []
        self.mapping        = {}
        self._plot_refs     = []
        self._event_markers = []
        self._build_ui()
        self._apply_dark_theme()

    def _build_ui(self):
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        root = QtWidgets.QHBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(8)

        
        # Right
        right_panel = QtWidgets.QVBoxLayout()
        right_panel.setSpacing(8)
        rl = QtWidgets.QLabel("3D Orientation")
        rl.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        rl.setStyleSheet("font-weight:bold; font-size:13px; color:#cdd6f4;")
        right_panel.addWidget(rl)
        self.rocket_view = RocketView()
        self.rocket_view.setMinimumSize(300, 320)
        right_panel.addWidget(self.rocket_view, stretch=3)
        rw = QtWidgets.QWidget()
        rw.setLayout(right_panel)
        root.addWidget(rw)

    def _apply_dark_theme(self):
        pg.setConfigOption("background", "#1e1e2e")
        pg.setConfigOption("foreground", "#cdd6f4")
        self.setStyleSheet(
            "QMainWindow, QWidget { background:#1e1e2e; color:#cdd6f4; }"
            "QListWidget { background:#181825; border:1px solid #45475a; border-radius:4px; font-size:12px; }"
            "QListWidget::item { padding:4px 8px; }"
            "QListWidget::item:hover { background:#313244; }"
            "QPushButton { background:#313244; border:1px solid #45475a; border-radius:4px; padding:4px 12px; color:#cdd6f4; }"
            "QPushButton:hover { background:#45475a; }"
            "QPushButton:checked { background:#45475a; border-color:#89b4fa; color:#89b4fa; }"
            "QLineEdit, QComboBox, QDoubleSpinBox { background:#181825; border:1px solid #45475a; border-radius:4px; padding:4px; color:#cdd6f4; }"
            "QComboBox::drop-down { border:none; }"
            "QGroupBox { border:1px solid #45475a; border-radius:4px; margin-top:8px; padding-top:6px; }"
            "QGroupBox::title { subcontrol-origin:margin; left:8px; color:#a6adc8; }"
            "QDialogButtonBox QPushButton { min-width:80px; }"
        )

    def _load_csv(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open Flight Log", "", "CSV Files (*.csv);;All Files (*)")
        if path:
            self._parse_csv(path)

    def _parse_csv(self, path):
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            headers = reader.fieldnames or []
            rows    = list(reader)
        if not rows or not headers:
            QtWidgets.QMessageBox.warning(self, "Error", "Empty or invalid CSV.")
            return
        self.data = {}
        for h in headers:
            try:
                self.data[h] = np.array([float(r[h]) if r[h] is not None and r[h] != '' else float('nan') for r in rows])
            except (ValueError, KeyError):
                pass
        self.headers = list(self.data.keys())
        if not self.headers:
            QtWidgets.QMessageBox.warning(self, "Error", "No numeric columns found.")
            return
        self.mapping = {}        # reset mapping, open dialog to set it
        self._raw_data = {k: v.copy() for k, v in self.data.items()}  # save original
        self._open_mapping()

    def _open_mapping(self):
        if not self.headers:
            QtWidgets.QMessageBox.information(self, "No Data", "Load a CSV file first.")
            return
        dlg = MappingDialog(self.headers, self.mapping, self)
        if dlg.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            self.mapping = dlg.get_mapping()
            self._trim_to_flight()
            self._refresh_all()

    def _refresh_all(self):
        time_col = self.mapping.get("time")
        if time_col:
            self.plot_widget.setLabel("bottom", time_col)
        self._populate_list()
        self._rebuild_plots()
        self._update_motor_panel()
        self._update_events_bar()

    def _populate_list(self):
        time_col = self.mapping.get("time")
        self.channel_list.blockSignals(True)
        self.channel_list.clear()
        for i, h in enumerate(self.headers):
            if h == time_col:
                continue
            item = QtWidgets.QListWidgetItem(h)
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(QtCore.Qt.CheckState.Unchecked)
            item.setForeground(QtGui.QColor(color_for_index(i)))
            self.channel_list.addItem(item)
        self.channel_list.blockSignals(False)

    def _on_item_changed(self, _): self._rebuild_plots()

    def _filter_list(self, text):
        for i in range(self.channel_list.count()):
            item = self.channel_list.item(i)
            item.setHidden(text.lower() not in item.text().lower())

    def _select_all(self):
        self.channel_list.blockSignals(True)
        for i in range(self.channel_list.count()):
            self.channel_list.item(i).setCheckState(QtCore.Qt.CheckState.Checked)
        self.channel_list.blockSignals(False)
        self._rebuild_plots()

    def _clear_all(self):
        self.channel_list.blockSignals(True)
        for i in range(self.channel_list.count()):
            self.channel_list.item(i).setCheckState(QtCore.Qt.CheckState.Unchecked)
        self.channel_list.blockSignals(False)
        self._rebuild_plots()

    def _checked_channels(self):
        checked = []
        for i in range(self.channel_list.count()):
            item = self.channel_list.item(i)
            if item.checkState() == QtCore.Qt.CheckState.Checked:
                checked.append((i, item.text()))
        return checked

    def _rebuild_plots(self):
        self.plot_widget.clear()
        self.legend = self.plot_widget.addLegend(offset=(10, 5))
        self.plot_widget.addItem(self.vline, ignoreBounds=True)
        self._plot_refs     = []
        self._event_markers = []
        time_col = self.mapping.get("time")
        if not time_col or time_col not in self.data:
            return
        t = self.data[time_col]
        for list_idx, ch in self._checked_channels():
            pen = pg.mkPen(color=color_for_index(list_idx), width=1.5)
            self._plot_refs.append(self.plot_widget.plot(t, self.data[ch], pen=pen, name=ch))
        if self.events_bar.display_btn.isChecked():
            self._on_markers_toggled(True, self.events_bar._events)

    def _reset_zoom(self):
        self.plot_widget.autoRange()

    def _update_motor_panel(self):
        time_col = self.mapping.get("time")
        az_col   = self.mapping.get("accel_world_z")
        mass_kg  = self.mapping.get("mass_kg", 1.0) or 1.0
        if time_col and az_col and time_col in self.data and az_col in self.data:
            accel  = self.data[az_col]
            thrust = derive_thrust(accel, mass_kg)
            stats  = compute_motor_stats(self.data[time_col], thrust, accel_world_z=accel)
            self.motor_panel.update_stats(stats)
        else:
            self.motor_panel.clear()

    def _update_events_bar(self):
        time_col = self.mapping.get("time")
        az_col   = self.mapping.get("accel_world_z")
        alt_col  = self.mapping.get("altitude")
        if not time_col or time_col not in self.data:
            self.events_bar.clear()
            return
        accel  = self.data[az_col]  if (az_col  and az_col  in self.data) else None
        alt    = self.data[alt_col] if (alt_col and alt_col in self.data) else None
        events = compute_flight_events(self.data[time_col], accel_world_z=accel, altitude=alt)
        self.events_bar.update_events(events)

    def _on_markers_toggled(self, show, events):
        for m in self._event_markers:
            self.plot_widget.removeItem(m)
        self._event_markers.clear()
        if not show:
            return
        time_col  = self.mapping.get("time")
        launch_ms = events.get("launch_ms")
        if not time_col or time_col not in self.data or launch_ms is None:
            return

        # (raw_ms_position, label, color)
        marker_defs = [
            (launch_ms,                                  "Launch",   "#e74c3c") if launch_ms                 is not None else None,
            (launch_ms + events["burn_time"] * 1000.0,  "Burnout",  "#f39c12") if events.get("burn_time")   is not None else None,
            (events.get("apogee_ms"),                    "Apogee",   "#2ecc71") if events.get("apogee_ms")   is not None else None,
            (events.get("land_ms"),                      "Landing",  "#3498db") if events.get("land_ms")     is not None else None,
        ]
        for entry in marker_defs:
            if entry is None:
                continue
            raw_t, label, color = entry
            if raw_t is None:
                continue
            line = pg.InfiniteLine(
                pos=raw_t, angle=90, movable=False,
                pen=pg.mkPen(color, width=1.5, style=QtCore.Qt.PenStyle.DashLine),
                label=label,
                labelOpts={"color": color, "position": 0.92,
                            "anchors": [(0.5, 0), (0.5, 1)]}
            )
            self.plot_widget.addItem(line)
            self._event_markers.append(line)

#Here is where the quaternions feed into the 3d model
    def _update_readout_and_3d(self, idx, t):
        parts = ["<b>t = %.4f</b>" % t[idx]]
        for list_idx, ch in self._checked_channels():
            val   = self.data[ch][idx]
            color = color_for_index(list_idx)
            parts.append('<span style="color:%s;">|%s: <b>%.4f</b></span>' % (color, ch, val))
        self.readout_label.setText("  ".join(parts))
        cols = [self.mapping.get(k) for k in ("qw", "qx", "qy", "qz")]
        if all(c and c in self.data for c in cols):
            self.rocket_view.set_quaternion(
                float(self.data[cols[0]][idx]), float(self.data[cols[1]][idx]),
                float(self.data[cols[2]][idx]), float(self.data[cols[3]][idx]),
            )

    def _on_mouse_moved(self, pos):
        vb = self.plot_widget.getViewBox()
        if not self.plot_widget.sceneBoundingRect().contains(pos):
            return
        x = vb.mapSceneToView(pos).x()
        self.vline.setPos(x)
        time_col = self.mapping.get("time")
        if not time_col or time_col not in self.data:
            return
        t   = self.data[time_col]
        idx = int(np.clip(np.searchsorted(t, x), 0, len(t) - 1))
        self._update_readout_and_3d(idx, t)
    def _trim_to_flight(self):
        time_col = self.mapping.get("time")
        az_col   = self.mapping.get("accel_world_z")
        if not time_col or time_col not in self.data:
            return

        time_ms   = self.data[time_col]
        accel     = self.data[az_col] if (az_col and az_col in self.data) else None
        alt_col   = self.mapping.get("altitude")
        altitude  = self.data[alt_col] if (alt_col and alt_col in self.data) else None

        events = compute_flight_events(time_ms, accel_world_z=accel, altitude=altitude)

        launch_ms = events.get("launch_ms")
        land_ms   = events.get("land_ms")

        if launch_ms is None:
            return  # can't trim without knowing launch

        t_start = launch_ms - 1000.0                                          # 1s before launch
        t_end   = (land_ms + 1000.0) if land_ms is not None else time_ms[-1] # 1s after landing

        mask = (time_ms >= t_start) & (time_ms <= t_end)
        if not np.any(mask):
            return

        for col in list(self.data.keys()):
            self.data[col] = self.data[col][mask]


def main():
    app = QtWidgets.QApplication(sys.argv)
    app.setApplicationName("Rocket Flight Reviewer")
    window = FlightReviewer()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()