import sys
import csv
import numpy as np
from PyQt6 import QtWidgets, QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl

from motor import derive_thrust, compute_motor_stats
from events_bar import FlightEventsBar, compute_flight_events

QtWidgets.QApplication.setAttribute(QtCore.Qt.ApplicationAttribute.AA_ShareOpenGLContexts)

fmt = QtGui.QSurfaceFormat()
fmt.setVersion(4, 1)
fmt.setProfile(QtGui.QSurfaceFormat.OpenGLContextProfile.CoreProfile)
fmt.setDepthBufferSize(24)
QtGui.QSurfaceFormat.setDefaultFormat(fmt)

THEME = {
    "bg": "#080C14",
    "panel": "#101724",
    "panel_2": "#131D2B",
    "panel_3": "#182235",
    "border": "#263245",
    "border_strong": "#34445C",
    "text": "#E6EDF7",
    "muted": "#91A0B5",
    "subtle": "#5F6F86",
    "primary": "#38BDF8",
    "success": "#22C55E",
    "warning": "#F59E0B",
    "danger": "#EF4444",
    "accent": "#A78BFA",
}

COLORS = [
    "#38BDF8", "#22C55E", "#F59E0B", "#A78BFA", "#F43F5E",
    "#14B8A6", "#EAB308", "#60A5FA", "#FB7185", "#34D399",
    "#F97316", "#818CF8", "#2DD4BF", "#C084FC", "#FACC15",
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
        "path_x":        ["path_x", "pos_x", "position_x", "x", "east", "easting", "longitude", "lon"],
        "path_y":        ["path_y", "pos_y", "position_y", "y", "north", "northing", "latitude", "lat"],
        "path_z":        ["path_z", "pos_z", "position_z", "z", "up", "altitude", "altitude_m", "height"],
    }
    LABELS = {
        "time":          "Time axis",
        "qw":            "Quaternion W",
        "qx":            "Quaternion X",
        "qy":            "Quaternion Y",
        "qz":            "Quaternion Z",
        "accel_world_z": "Accel World Z (m/s^2)",
        "altitude":      "Altitude (m)",
        "path_x":        "Path X / East (optional)",
        "path_y":        "Path Y / North (optional)",
        "path_z":        "Path Z / Up (optional)",
    }
    GROUPS = [
        ("General",    ["time"]),
        ("Quaternion", ["qw", "qx", "qy", "qz"]),
        ("Motor",      ["accel_world_z"]),
        ("3D Flight",  ["altitude", "path_x", "path_y", "path_z"]),
    ]

    def __init__(self, headers, current_mapping, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Column Mapping")
        self.setMinimumWidth(520)
        self.headers   = headers
        self.combos    = {}
        self.mass_spin = None
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(14)
        title = QtWidgets.QLabel("Column Mapping")
        title.setObjectName("DialogTitle")
        layout.addWidget(title)
        hint = QtWidgets.QLabel("Map CSV columns to program fields. Leave '-- none --' to skip.")
        hint.setObjectName("MutedText")
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
        self.setObjectName("RocketViewport")
        self.setCameraPosition(distance=5, elevation=20, azimuth=45)
        self.setBackgroundColor(QtGui.QColor(THEME["bg"]))
        grid = gl.GLGridItem()
        grid.setSize(6, 6)
        grid.setSpacing(1, 1)
        grid.setColor((52, 68, 92, 120))
        self.addItem(grid)
        self.addItem(gl.GLAxisItem(size=QtGui.QVector3D(1.0, 1.0, 2.5)))
        verts, faces, face_colors = build_rocket_mesh()
        md = gl.MeshData(vertexes=verts, faces=faces, faceColors=face_colors)
        self.rocket = gl.GLMeshItem(meshdata=md, smooth=True, shader="normalColor", glOptions="opaque")
        self.addItem(self.rocket)
        self.info_label = QtWidgets.QLabel("qw=1.000  qx=0.000  qy=0.000  qz=0.000", self)
        self.info_label.setStyleSheet(
            "color:#E6EDF7; background:rgba(16,23,36,220); "
            "font-family:'Menlo','Monaco',monospace; font-size:11px; "
            "padding:6px 10px; border:1px solid rgba(56,189,248,80); border-radius:6px;"
        )
        self.info_label.move(10, 10)
        self.info_label.adjustSize()

    def set_quaternion(self, w, x, y, z):
        q = QtGui.QQuaternion(w, x, y, z).normalized()
        axis, angle = q.getAxisAndAngle()
        self.rocket.resetTransform()
        self.rocket.rotate(angle, axis.x(), axis.y(), axis.z())
        self.info_label.setText("qw=%+.3f  qx=%+.3f  qy=%+.3f  qz=%+.3f" % (w, x, y, z))
        self.info_label.adjustSize()


class FlightPathView(gl.GLViewWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("RocketViewport")
        self.setCameraPosition(distance=9, elevation=24, azimuth=45)
        self.setBackgroundColor(QtGui.QColor(THEME["bg"]))

        grid = gl.GLGridItem()
        grid.setSize(8, 8)
        grid.setSpacing(1, 1)
        grid.setColor((52, 68, 92, 120))
        self.addItem(grid)
        self.addItem(gl.GLAxisItem(size=QtGui.QVector3D(1.2, 1.2, 2.0)))

        self.path_line = gl.GLLinePlotItem(
            pos=np.zeros((0, 3), dtype=np.float32),
            color=(0.22, 0.74, 0.97, 1.0),
            width=2.0,
            antialias=True,
            mode="line_strip",
        )
        self.addItem(self.path_line)
        self.path_marker = gl.GLScatterPlotItem(
            pos=np.zeros((0, 3), dtype=np.float32),
            color=(0.96, 0.62, 0.04, 1.0),
            size=10.0,
            pxMode=True,
        )
        self.addItem(self.path_marker)

        self.path_label = QtWidgets.QLabel("Path unavailable: map X/Y, Z, or altitude", self)
        self.path_label.setStyleSheet(
            "color:#91A0B5; background:rgba(16,23,36,220); "
            "font-size:11px; padding:6px 10px; "
            "border:1px solid rgba(52,68,92,170); border-radius:6px;"
        )
        self.path_label.move(10, 10)
        self.path_label.adjustSize()

    def set_path_status(self, text):
        self.path_label.setText(text)
        self.path_label.adjustSize()

    def clear_path(self):
        self.path_line.setData(pos=np.zeros((0, 3), dtype=np.float32))
        self.path_marker.setData(pos=np.zeros((0, 3), dtype=np.float32))
        self.set_path_status("Path unavailable: map X/Y, Z, or altitude")

    def set_flight_path(self, points, idx):
        if points is None or len(points) == 0:
            self.clear_path()
            return
        idx = int(np.clip(idx, 0, len(points) - 1))
        trail = np.asarray(points[:idx + 1], dtype=np.float32)
        current = np.asarray(points[idx], dtype=np.float32)
        self.path_line.setData(pos=trail)
        self.path_marker.setData(pos=current.reshape(1, 3))


class MotorPanel(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("Panel")
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(14, 14, 14, 14)
        layout.setSpacing(10)
        title = QtWidgets.QLabel("Motor Classification")
        title.setObjectName("SectionTitle")
        layout.addWidget(title)
        self.designation_label = QtWidgets.QLabel("--")
        self.designation_label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.designation_label.setStyleSheet(
            "font-size:46px; font-weight:800; color:#F59E0B; "
            "background:#080C14; border:1px solid #34445C; border-radius:8px; padding:12px;"
        )
        layout.addWidget(self.designation_label)
        grid = QtWidgets.QGridLayout()
        grid.setHorizontalSpacing(10)
        grid.setVerticalSpacing(8)
        self.stat_labels = {}
        for row, (key, label, unit) in enumerate([
            ("total_impulse", "Total Impulse", "N*s"),
            ("burn_time",     "Burn Time",     "s"),
            ("avg_thrust",    "Avg Thrust",    "N"),
            ("peak_thrust",   "Peak Thrust",   "N"),
        ]):
            nl = QtWidgets.QLabel(label)
            nl.setObjectName("MetricLabel")
            vl = QtWidgets.QLabel("--")
            vl.setStyleSheet(
                "color:#E6EDF7; font-family:'Menlo','Monaco',monospace; "
                "font-size:12px; font-weight:700; background:#080C14; "
                "border:1px solid #263245; border-radius:5px; padding:4px 8px;"
            )
            ul = QtWidgets.QLabel(unit)
            ul.setObjectName("UnitLabel")
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
        self.resize(1760, 940)
        self.data           = {}
        self.headers        = []
        self.mapping        = {}
        self._plot_refs     = []
        self._event_markers = []
        self._measure_times = []
        self._measure_markers = []
        self._path_points = None
        self._path_source = None
        self._build_ui()
        self._apply_dark_theme()

    def _panel(self, name=None):
        panel = QtWidgets.QWidget()
        panel.setObjectName(name or "Panel")
        return panel

    def _section_title(self, text, detail=None):
        row = QtWidgets.QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        title = QtWidgets.QLabel(text)
        title.setObjectName("SectionTitle")
        row.addWidget(title)
        if detail:
            row.addStretch()
            detail_label = QtWidgets.QLabel(detail)
            detail_label.setObjectName("SectionDetail")
            row.addWidget(detail_label)
        return row

    def _build_ui(self):
        central = QtWidgets.QWidget()
        central.setObjectName("AppRoot")
        self.setCentralWidget(central)
        root = QtWidgets.QHBoxLayout(central)
        root.setContentsMargins(14, 14, 14, 14)
        root.setSpacing(12)

        # Left
        left_panel = self._panel()
        left = QtWidgets.QVBoxLayout(left_panel)
        left.setContentsMargins(14, 14, 14, 14)
        left.setSpacing(12)
        left.addLayout(self._section_title("Flight Data", "CSV"))
        self.loaded_file_label = QtWidgets.QLabel("No file loaded")
        self.loaded_file_label.setObjectName("StatusPill")
        self.loaded_file_label.setToolTip("Current flight log")
        self.loaded_file_label.setWordWrap(True)
        left.addWidget(self.loaded_file_label)

        load_btn = QtWidgets.QPushButton("Load CSV")
        load_btn.setObjectName("PrimaryButton")
        load_btn.setFixedHeight(38)
        load_btn.clicked.connect(self._load_csv)
        left.addWidget(load_btn)

        map_btn = QtWidgets.QPushButton("Column Mapping")
        map_btn.setFixedHeight(34)
        map_btn.clicked.connect(self._open_mapping)
        left.addWidget(map_btn)

        left.addSpacing(2)
        left.addLayout(self._section_title("Channels"))
        self.search_bar = QtWidgets.QLineEdit()
        self.search_bar.setPlaceholderText("Filter channels...")
        self.search_bar.setClearButtonEnabled(True)
        self.search_bar.textChanged.connect(self._filter_list)
        left.addWidget(self.search_bar)

        self.channel_list = QtWidgets.QListWidget()
        self.channel_list.setObjectName("ChannelList")
        self.channel_list.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.NoSelection)
        self.channel_list.itemChanged.connect(self._on_item_changed)
        left.addWidget(self.channel_list, stretch=1)

        self.selection_status = QtWidgets.QLabel("0 selected")
        self.selection_status.setObjectName("MutedText")
        left.addWidget(self.selection_status)

        btn_row = QtWidgets.QHBoxLayout()
        btn_row.setSpacing(8)
        for label, slot in [("Select All", self._select_all), ("Clear All", self._clear_all)]:
            b = QtWidgets.QPushButton(label)
            b.setFixedHeight(32)
            b.clicked.connect(slot)
            btn_row.addWidget(b)
        left.addLayout(btn_row)
        left_panel.setFixedWidth(270)
        root.addWidget(left_panel)

        # Center
        center_panel = self._panel()
        center = QtWidgets.QVBoxLayout(center_panel)
        center.setContentsMargins(14, 14, 14, 14)
        center.setSpacing(10)
        tb = QtWidgets.QHBoxLayout()
        tb.setSpacing(10)
        title_stack = QtWidgets.QVBoxLayout()
        title_stack.setSpacing(2)
        chart_title = QtWidgets.QLabel("Telemetry Plot")
        chart_title.setObjectName("SectionTitle")
        chart_subtitle = QtWidgets.QLabel("Hover the graph to inspect synchronized channel values and 3D attitude.")
        chart_subtitle.setObjectName("MutedText")
        title_stack.addWidget(chart_title)
        title_stack.addWidget(chart_subtitle)
        tb.addLayout(title_stack)
        rb = QtWidgets.QPushButton("Reset Zoom")
        rb.setFixedHeight(32)
        rb.clicked.connect(self._reset_zoom)
        tb.addStretch()
        tb.addWidget(rb)
        center.addLayout(tb)

        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setObjectName("PlotWidget")
        self.plot_widget.showGrid(x=True, y=True, alpha=0.3)
        self.plot_widget.setLabel("bottom", "Time")
        self.plot_widget.setMenuEnabled(False)
        self.legend = self.plot_widget.addLegend(offset=(10, 5))
        center.addWidget(self.plot_widget, stretch=1)

        self.readout_label = QtWidgets.QLabel("")
        self.readout_label.setWordWrap(True)
        self.readout_label.setObjectName("Readout")
        self.readout_label.setText("Load a CSV, map columns, then select channels to inspect telemetry.")
        self.readout_label.setContentsMargins(10, 7, 10, 7)
        self.readout_label.setStyleSheet(
            "font-family:'Menlo','Monaco',monospace; font-size:11px; "
            "background:#080C14; border:1px solid #263245; border-radius:7px; padding:6px;"
        )
        self.readout_label.setFixedHeight(58)
        center.addWidget(self.readout_label)

        measure_row = QtWidgets.QHBoxLayout()
        measure_row.setSpacing(8)
        self.measure_label = QtWidgets.QLabel("Click the plot twice to measure time difference.")
        self.measure_label.setObjectName("StatusPill")
        self.measure_label.setWordWrap(True)
        measure_row.addWidget(self.measure_label, stretch=1)
        clear_measure_btn = QtWidgets.QPushButton("Clear Measure")
        clear_measure_btn.setFixedHeight(32)
        clear_measure_btn.clicked.connect(self._clear_measurement)
        measure_row.addWidget(clear_measure_btn)
        center.addLayout(measure_row)

        self.events_bar = FlightEventsBar()
        self.events_bar.markers_toggled.connect(self._on_markers_toggled)
        center.addWidget(self.events_bar)
        root.addWidget(center_panel, stretch=3)

        # Right
        right_shell = self._panel()
        right_panel = QtWidgets.QVBoxLayout(right_shell)
        right_panel.setContentsMargins(14, 14, 14, 14)
        right_panel.setSpacing(12)
        right_panel.addLayout(self._section_title("3D Views", "Replay"))
        self.visual_tabs = QtWidgets.QTabWidget()
        self.visual_tabs.setObjectName("VisualTabs")
        self.rocket_view = RocketView()
        self.flight_path_view = FlightPathView()
        self.rocket_view.setMinimumSize(320, 320)
        self.flight_path_view.setMinimumSize(320, 320)
        self.visual_tabs.addTab(self.rocket_view, "Orientation")
        self.visual_tabs.addTab(self.flight_path_view, "Flight Path")
        right_panel.addWidget(self.visual_tabs, stretch=3)
        self.motor_panel = MotorPanel()
        self.motor_panel.setMinimumHeight(230)
        right_panel.addWidget(self.motor_panel, stretch=2)
        right_shell.setFixedWidth(360)
        root.addWidget(right_shell)

        self.vline = pg.InfiniteLine(angle=90, movable=False,
            pen=pg.mkPen(THEME["text"], width=1, style=QtCore.Qt.PenStyle.DashLine))
        self.plot_widget.addItem(self.vline, ignoreBounds=True)
        self.plot_widget.scene().sigMouseMoved.connect(self._on_mouse_moved)
        self.plot_widget.scene().sigMouseClicked.connect(self._on_plot_clicked)

    def _apply_dark_theme(self):
        pg.setConfigOption("background", THEME["panel"])
        pg.setConfigOption("foreground", THEME["text"])
        self.plot_widget.setBackground(THEME["panel"])
        self.plot_widget.getAxis("bottom").setPen(pg.mkPen(THEME["border_strong"]))
        self.plot_widget.getAxis("left").setPen(pg.mkPen(THEME["border_strong"]))
        self.plot_widget.getAxis("bottom").setTextPen(pg.mkPen(THEME["muted"]))
        self.plot_widget.getAxis("left").setTextPen(pg.mkPen(THEME["muted"]))
        self.setStyleSheet(
            f"""
            QMainWindow, QWidget#AppRoot {{
                background:{THEME["bg"]};
                color:{THEME["text"]};
                font-family:-apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
                font-size:12px;
            }}
            QWidget#Panel {{
                background:{THEME["panel"]};
                border:1px solid {THEME["border"]};
                border-radius:8px;
            }}
            QWidget#RocketViewport {{
                border:1px solid {THEME["border"]};
                border-radius:8px;
            }}
            QTabWidget#VisualTabs::pane {{
                border:1px solid {THEME["border"]};
                border-radius:8px;
                top:-1px;
            }}
            QTabWidget#VisualTabs QTabBar::tab {{
                background:{THEME["bg"]};
                border:1px solid {THEME["border"]};
                border-bottom:none;
                border-top-left-radius:6px;
                border-top-right-radius:6px;
                color:{THEME["muted"]};
                padding:7px 12px;
                font-weight:650;
            }}
            QTabWidget#VisualTabs QTabBar::tab:selected {{
                background:{THEME["panel_2"]};
                color:{THEME["text"]};
                border-color:{THEME["border_strong"]};
            }}
            QTabWidget#VisualTabs QTabBar::tab:hover {{
                color:{THEME["primary"]};
            }}
            QLabel#SectionTitle {{
                background:transparent;
                color:{THEME["text"]};
                font-size:13px;
                font-weight:700;
                letter-spacing:0px;
            }}
            QLabel#SectionDetail, QLabel#MutedText, QLabel#MetricLabel {{
                background:transparent;
                color:{THEME["muted"]};
                font-size:11px;
            }}
            QLabel#UnitLabel {{
                background:transparent;
                color:{THEME["subtle"]};
                font-size:11px;
            }}
            QLabel#DialogTitle {{
                background:transparent;
                color:{THEME["text"]};
                font-size:18px;
                font-weight:800;
            }}
            QLabel#StatusPill {{
                background:{THEME["panel_2"]};
                border:1px solid {THEME["border"]};
                border-radius:7px;
                color:{THEME["muted"]};
                padding:8px 10px;
                font-family:"Menlo", "Monaco", monospace;
                font-size:11px;
            }}
            QLabel#Readout {{
                color:{THEME["text"]};
            }}
            QListWidget#ChannelList {{
                background:{THEME["bg"]};
                border:1px solid {THEME["border"]};
                border-radius:7px;
                font-size:12px;
                outline:0;
            }}
            QListWidget#ChannelList::item {{
                padding:6px 8px;
                border-bottom:1px solid rgba(38,50,69,90);
            }}
            QListWidget#ChannelList::item:hover {{
                background:{THEME["panel_3"]};
            }}
            QListWidget#ChannelList::indicator {{
                width:14px;
                height:14px;
            }}
            QPushButton {{
                background:{THEME["panel_2"]};
                border:1px solid {THEME["border_strong"]};
                border-radius:6px;
                padding:5px 12px;
                color:{THEME["text"]};
                font-weight:650;
            }}
            QPushButton:hover {{
                background:{THEME["panel_3"]};
                border-color:{THEME["primary"]};
            }}
            QPushButton:pressed {{
                background:{THEME["bg"]};
            }}
            QPushButton:focus {{
                border:1px solid {THEME["primary"]};
            }}
            QPushButton#PrimaryButton {{
                background:{THEME["primary"]};
                color:#04111E;
                border-color:{THEME["primary"]};
            }}
            QPushButton#PrimaryButton:hover {{
                background:#7DD3FC;
                border-color:#7DD3FC;
            }}
            QPushButton:checked {{
                background:#0E7490;
                border-color:{THEME["primary"]};
                color:#ECFEFF;
            }}
            QLineEdit, QComboBox, QDoubleSpinBox {{
                background:{THEME["bg"]};
                border:1px solid {THEME["border"]};
                border-radius:6px;
                padding:6px 8px;
                color:{THEME["text"]};
                selection-background-color:{THEME["primary"]};
                selection-color:#04111E;
            }}
            QLineEdit:focus, QComboBox:focus, QDoubleSpinBox:focus {{
                border:1px solid {THEME["primary"]};
            }}
            QComboBox::drop-down {{
                border:none;
                width:24px;
            }}
            QGroupBox {{
                border:1px solid {THEME["border"]};
                border-radius:7px;
                margin-top:10px;
                padding:10px 8px 8px 8px;
                color:{THEME["text"]};
                font-weight:700;
            }}
            QGroupBox::title {{
                subcontrol-origin:margin;
                left:10px;
                padding:0 5px;
                color:{THEME["muted"]};
            }}
            QDialog {{
                background:{THEME["panel"]};
                color:{THEME["text"]};
            }}
            QDialogButtonBox QPushButton {{
                min-width:88px;
            }}
            QToolTip {{
                background:{THEME["panel_3"]};
                color:{THEME["text"]};
                border:1px solid {THEME["border_strong"]};
                padding:6px;
            }}
            """
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
        file_name = QtCore.QFileInfo(path).fileName()
        self.loaded_file_label.setText("%s\n%d ch  |  %d rows" % (
            file_name, len(self.headers), len(rows)
        ))
        self.loaded_file_label.setToolTip(path)
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
        self._update_flight_path()
        self._rebuild_plots()
        self._update_motor_panel()
        self._update_events_bar()

    def _normalize_path_points(self, x, y, z):
        points = np.column_stack([x, y, z]).astype(float)
        finite = np.all(np.isfinite(points), axis=1)
        if not np.any(finite):
            return None
        cleaned = points.copy()
        for axis in range(3):
            values = cleaned[:, axis]
            fill = np.nanmedian(values[finite])
            values[~np.isfinite(values)] = fill
            cleaned[:, axis] = values
        mins = np.nanmin(cleaned[finite], axis=0)
        maxs = np.nanmax(cleaned[finite], axis=0)
        cleaned -= (mins + maxs) / 2.0
        span = np.nanmax(maxs - mins)
        if not np.isfinite(span) or span <= 0:
            span = 1.0
        return (cleaned / span * 4.2).astype(np.float32)

    def _derived_path_from_attitude(self, altitude):
        cols = [self.mapping.get(k) for k in ("qw", "qx", "qy", "qz")]
        if not all(c and c in self.data for c in cols):
            x = np.zeros_like(altitude, dtype=float)
            y = np.zeros_like(altitude, dtype=float)
            return x, y, "Altitude-only trail"

        qw, qx, qy, qz = (self.data[c].astype(float) for c in cols)
        forward_x = 2.0 * (qx * qz + qw * qy)
        forward_y = 2.0 * (qy * qz - qw * qx)
        dz = np.diff(altitude, prepend=altitude[0])
        step = np.maximum(np.abs(dz), 0.05)
        x = np.cumsum(forward_x * step)
        y = np.cumsum(forward_y * step)
        return x, y, "Derived trail from altitude + attitude"

    def _update_flight_path(self):
        x_col = self.mapping.get("path_x")
        y_col = self.mapping.get("path_y")
        z_col = self.mapping.get("path_z")
        alt_col = self.mapping.get("altitude")

        if x_col and y_col and x_col in self.data and y_col in self.data:
            x = self.data[x_col].astype(float)
            y = self.data[y_col].astype(float)
            if z_col and z_col in self.data:
                z = self.data[z_col].astype(float)
                source = "Mapped X/Y/Z trail"
            elif alt_col and alt_col in self.data:
                z = self.data[alt_col].astype(float)
                source = "Mapped X/Y + altitude trail"
            else:
                z = np.zeros_like(x, dtype=float)
                source = "Mapped X/Y flat trail"
            self._path_points = self._normalize_path_points(x, y, z)
        else:
            if not alt_col or alt_col not in self.data:
                self._path_points = None
                self._path_source = None
                self.flight_path_view.clear_path()
                return
            altitude = self.data[alt_col].astype(float)
            x, y, source = self._derived_path_from_attitude(altitude)
            self._path_points = self._normalize_path_points(x, y, altitude)

        self._path_source = source if self._path_points is not None else None
        if self._path_points is None:
            self.flight_path_view.clear_path()
        else:
            self.flight_path_view.set_path_status("%s: hover plot to replay path" % source)

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
            item.setToolTip(h)
            item.setData(QtCore.Qt.ItemDataRole.UserRole, i)
            item.setForeground(QtGui.QColor(color_for_index(i)))
            self.channel_list.addItem(item)
        self.channel_list.blockSignals(False)
        self._update_selection_status()

    def _on_item_changed(self, _):
        self._update_selection_status()
        self._rebuild_plots()

    def _filter_list(self, text):
        for i in range(self.channel_list.count()):
            item = self.channel_list.item(i)
            item.setHidden(text.lower() not in item.text().lower())
        self._update_selection_status()

    def _select_all(self):
        self.channel_list.blockSignals(True)
        for i in range(self.channel_list.count()):
            self.channel_list.item(i).setCheckState(QtCore.Qt.CheckState.Checked)
        self.channel_list.blockSignals(False)
        self._update_selection_status()
        self._rebuild_plots()

    def _clear_all(self):
        self.channel_list.blockSignals(True)
        for i in range(self.channel_list.count()):
            self.channel_list.item(i).setCheckState(QtCore.Qt.CheckState.Unchecked)
        self.channel_list.blockSignals(False)
        self._update_selection_status()
        self._rebuild_plots()

    def _checked_channels(self):
        checked = []
        for i in range(self.channel_list.count()):
            item = self.channel_list.item(i)
            if item.checkState() == QtCore.Qt.CheckState.Checked:
                color_idx = item.data(QtCore.Qt.ItemDataRole.UserRole)
                checked.append((int(color_idx), item.text()))
        return checked

    def _update_selection_status(self):
        total = self.channel_list.count()
        selected = len(self._checked_channels())
        visible = sum(
            1 for i in range(total)
            if not self.channel_list.item(i).isHidden()
        )
        self.selection_status.setText("%d selected  |  %d visible  |  %d total" % (
            selected, visible, total
        ))

    def _rebuild_plots(self):
        self.plot_widget.clear()
        self.legend = self.plot_widget.addLegend(offset=(10, 5))
        self.plot_widget.addItem(self.vline, ignoreBounds=True)
        self._plot_refs     = []
        self._event_markers = []
        self._measure_markers = []
        time_col = self.mapping.get("time")
        if not time_col or time_col not in self.data:
            return
        t = self.data[time_col]
        for list_idx, ch in self._checked_channels():
            pen = pg.mkPen(color=color_for_index(list_idx), width=1.5)
            self._plot_refs.append(self.plot_widget.plot(t, self.data[ch], pen=pen, name=ch))
        if self.events_bar.display_btn.isChecked():
            self._on_markers_toggled(True, self.events_bar._events)
        self._redraw_measurement_markers()

    def _reset_zoom(self):
        self.plot_widget.autoRange()

    def _format_time_delta(self, dt):
        unit = self.mapping.get("time") or "time"
        if "ms" in unit.lower():
            return "%.3f ms  |  %.4f s" % (dt, dt / 1000.0)
        return "%.4f %s" % (dt, unit)

    def _clear_measurement(self):
        self._measure_times.clear()
        for marker in self._measure_markers:
            self.plot_widget.removeItem(marker)
        self._measure_markers.clear()
        self.measure_label.setText("Click the plot twice to measure time difference.")

    def _redraw_measurement_markers(self):
        for marker in self._measure_markers:
            self.plot_widget.removeItem(marker)
        self._measure_markers.clear()
        for idx, t_value in enumerate(self._measure_times):
            label = "A" if idx == 0 else "B"
            line = pg.InfiniteLine(
                pos=t_value,
                angle=90,
                movable=False,
                pen=pg.mkPen(THEME["accent"], width=2, style=QtCore.Qt.PenStyle.DotLine),
                label=label,
                labelOpts={
                    "color": THEME["accent"],
                    "position": 0.08,
                    "anchors": [(0.5, 0), (0.5, 1)],
                },
            )
            self.plot_widget.addItem(line)
            self._measure_markers.append(line)

    def _on_plot_clicked(self, event):
        if event.button() != QtCore.Qt.MouseButton.LeftButton:
            return
        if not self.plot_widget.sceneBoundingRect().contains(event.scenePos()):
            return
        time_col = self.mapping.get("time")
        if not time_col or time_col not in self.data:
            self.measure_label.setText("Load data and map a time column before measuring.")
            return

        vb = self.plot_widget.getViewBox()
        x = vb.mapSceneToView(event.scenePos()).x()
        t = self.data[time_col]
        idx = int(np.clip(np.searchsorted(t, x), 0, len(t) - 1))
        clicked_t = float(t[idx])

        if len(self._measure_times) >= 2:
            self._measure_times.clear()
        self._measure_times.append(clicked_t)
        self._redraw_measurement_markers()

        if len(self._measure_times) == 1:
            self.measure_label.setText("A = %.4f. Click a second point for delta." % clicked_t)
        else:
            a, b = self._measure_times
            dt = abs(b - a)
            self.measure_label.setText(
                "A = %.4f  |  B = %.4f  |  Δt = %s" % (
                    a, b, self._format_time_delta(dt)
                )
            )
        event.accept()

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

    def _update_readout_and_3d(self, idx, t):
        parts = ["<b>t = %.4f</b>" % t[idx]]
        for list_idx, ch in self._checked_channels():
            val   = self.data[ch][idx]
            color = color_for_index(list_idx)
            parts.append('<span style="color:%s;">|%s: <b>%.4f</b></span>' % (color, ch, val))
        self.readout_label.setText("  ".join(parts))
        if self._path_points is not None:
            self.flight_path_view.set_flight_path(self._path_points, idx)
            if self._path_source:
                self.flight_path_view.set_path_status("%s: point %d / %d" % (
                    self._path_source, idx + 1, len(self._path_points)
                ))
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
