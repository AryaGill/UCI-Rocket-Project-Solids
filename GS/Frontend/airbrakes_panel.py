from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QLabel, QFrame


class AirbrakesPanel(QWidget):
    ...

    def __init__(self, parent=None):
        super().__init__(parent)

        root = QVBoxLayout(self)
        root.setContentsMargins(6, 4, 0, 4)
        root.setSpacing(4)

        title = QLabel("Command Echo")
        title.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        title.setStyleSheet(
            "color: #b0b0b0; font-size: 10px; font-weight: bold;"
        )
        root.addWidget(title)

        # AB_Deployment value
        self.value_lbl = QLabel("AB_Deployment: --")
        self.value_lbl.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self.value_lbl.setStyleSheet(
            "color: #e0e0e0; font-size: 12px; font-family: monospace;"
        )
        root.addWidget(self.value_lbl)

        # NEW: predicted apogee
        self.pred_apo_lbl = QLabel("pred_apo: --")
        self.pred_apo_lbl.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self.pred_apo_lbl.setStyleSheet(
            "color: #e0e0e0; font-size: 12px; font-family: monospace;"
        )
        root.addWidget(self.pred_apo_lbl)

        # Command echo readout
        self.cmd_echo_lbl = QLabel("echo: --")
        self.cmd_echo_lbl.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self.cmd_echo_lbl.setStyleSheet(
            "color: #a855f7; font-size: 10px; font-family: monospace;"
        )
        root.addWidget(self.cmd_echo_lbl)

        self.status_lbl = QLabel("No Data")
        self.status_lbl.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self.status_lbl.setStyleSheet(
            "color: #808080; font-size: 10px;"
        )
        root.addWidget(self.status_lbl)

        self.setStyleSheet("background-color: #232323; border-radius: 6px;")
        self.setFixedHeight(130)   # taller to fit command echo line
        self.setMaximumWidth(260)

    def update_data(self, data: dict):
        """
        Expects telemetry dict with keys 'AB_Deployment', 'pred_apo', and 'command_echo'.
        """
        ab_value   = data.get("AB_Deployment", None)
        pred_apo   = data.get("pred_apo", None)
        cmd_echo   = data.get("command_echo", None)
    
        if ab_value is None and pred_apo is None:
            self.value_lbl.setText("AB_Deployment: --")
            self.pred_apo_lbl.setText("pred_apo: --")
            self.status_lbl.setText("No Data")
            return

        if ab_value is None:
            self.value_lbl.setText("AB_Deployment: --")
        else:
            self.value_lbl.setText(f"AB_Deployment: {ab_value}")

        if pred_apo is None:
            self.pred_apo_lbl.setText("pred_apo: --")
        else:
            self.pred_apo_lbl.setText(f"pred_apo: {pred_apo}")

        if cmd_echo is not None:
            self.cmd_echo_lbl.setText(f"echo: {cmd_echo}")

        self.status_lbl.setText("OK")