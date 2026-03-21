from __future__ import annotations

from PyQt6.QtWidgets import (QWidget, QHBoxLayout, QVBoxLayout,
                             QLabel, QFrame)
from PyQt6.QtCore import Qt


# ── Flight states that trigger each event ────────────────────────────────────
# "deploying" state  →  "deployed" state
_DROGUE_P_DEPLOYING = 4
_DROGUE_P_DEPLOYED  = 5
_DROGUE_S_DEPLOYING = 6
_DROGUE_S_DEPLOYED  = 7
_MAIN_P_DEPLOYING   = 8
_MAIN_P_DEPLOYED    = 9
_MAIN_S_DEPLOYING   = 10
_MAIN_S_DEPLOYED    = 11


class _DeployIndicator(QWidget):
    """
    Single parachute channel indicator.

    Three visual states
    ──────────────────
    idle      – grey,  "—"
    deploying – yellow, "DEPLOYING"
    deployed  – green,  "✓ DEPLOYED  @{time:.1f}s"
    """

    def __init__(self, label: str, accent: str, parent=None):
        super().__init__(parent)
        self._accent = accent

        root = QVBoxLayout(self)
        root.setContentsMargins(10, 6, 10, 6)
        root.setSpacing(3)

        # Channel name
        name = QLabel(label)
        name.setAlignment(Qt.AlignmentFlag.AlignCenter)
        name.setStyleSheet(
            f"color: {accent}; font-weight: bold; font-size: 10px;"
        )
        root.addWidget(name)

        # Big status text
        self.status_lbl = QLabel("—")
        self.status_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.status_lbl.setStyleSheet(
            "color: #606060; font-size: 13px; font-weight: bold; font-family: monospace;"
        )
        root.addWidget(self.status_lbl)

        # Timestamp / sub-text
        self.time_lbl = QLabel("")
        self.time_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.time_lbl.setStyleSheet(
            "color: #808080; font-size: 9px; font-family: monospace;"
        )
        root.addWidget(self.time_lbl)

        # Coloured bottom bar (indicator strip)
        self.bar = QFrame()
        self.bar.setFixedHeight(4)
        self.bar.setStyleSheet("background-color: #333333; border-radius: 2px;")
        root.addWidget(self.bar)

        self.setStyleSheet("background-color: #2a2a2a; border-radius: 5px;")
        self.setFixedWidth(130)

    # ── state setters ─────────────────────────────────────────────────────────

    def set_idle(self):
        self.status_lbl.setText("—")
        self.status_lbl.setStyleSheet(
            "color: #606060; font-size: 13px; font-weight: bold; font-family: monospace;"
        )
        self.time_lbl.setText("")
        self.bar.setStyleSheet("background-color: #333333; border-radius: 2px;")

    def set_deploying(self):
        self.status_lbl.setText("DEPLOYING")
        self.status_lbl.setStyleSheet(
            "color: #fffb00; font-size: 11px; font-weight: bold; font-family: monospace;"
        )
        self.time_lbl.setText("")
        self.bar.setStyleSheet("background-color: #fffb00; border-radius: 2px;")

    def set_deployed(self, timestamp: float | None = None):
        self.status_lbl.setText("✓  DEPLOYED")
        self.status_lbl.setStyleSheet(
            "color: #28dc5e; font-size: 11px; font-weight: bold; font-family: monospace;"
        )
        if timestamp is not None:
            self.time_lbl.setText(f"@ T+{timestamp:.1f} s")
        else:
            self.time_lbl.setText("")
        self.bar.setStyleSheet("background-color: #28dc5e; border-radius: 2px;")


class DeploymentStatusPanel(QWidget):
    """
    Horizontal panel showing deployment confirmation for all four parachute
    channels.  Feed it every telemetry packet via ``update_data(data)``.

    Channels tracked
    ────────────────
    Drogue Primary   – flight states 4 (deploying) → 5 (deployed)
    Drogue Secondary – flight states 6 (deploying) → 7 (deployed)
    Main Primary     – flight states 8 (deploying) → 9 (deployed)
    Main Secondary   – flight states 10 (deploying) → 11 (deployed)
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        self._prev_state: int | None = None

        outer = QVBoxLayout(self)
        outer.setContentsMargins(6, 4, 6, 4)
        outer.setSpacing(4)

        # Section title
        title = QLabel("Deployment Status")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title.setStyleSheet(
            "color: #b0b0b0; font-size: 10px; font-weight: bold;"
        )
        outer.addWidget(title)

        # Card row
        row = QHBoxLayout()
        row.setSpacing(8)
        row.setContentsMargins(0, 0, 0, 0)

        self._drogue_p = _DeployIndicator("Drogue\nPrimary",   "#ff6b35")
        self._drogue_s = _DeployIndicator("Drogue\nSecondary", "#ff9966")
        self._main_p   = _DeployIndicator("Main\nPrimary",     "#28dc5e")
        self._main_s   = _DeployIndicator("Main\nSecondary",   "#66ee99")

        for card in (self._drogue_p, self._drogue_s, self._main_p, self._main_s):
            row.addWidget(card)

        outer.addLayout(row)

        self.setStyleSheet("background-color: #232323; border-radius: 6px;")
        self.setFixedHeight(105)

    # ── public API ────────────────────────────────────────────────────────────

    def reset(self):
        """Reset all indicators to idle (call on graph clear / new flight)."""
        for card in (self._drogue_p, self._drogue_s, self._main_p, self._main_s):
            card.set_idle()
        self._prev_state = None

    def update_data(self, data: dict):
        """
        Feed a telemetry packet.  Reads 'flight_state' and 'Time'.
        Transitions are edge-triggered so indicators latch once deployed.
        """
        raw_state = data.get("flight_state")
        if raw_state is None:
            return

        try:
            state = int(raw_state)
        except (ValueError, TypeError):
            return

        timestamp = data.get("Time")  # may be None

        # Only act on state *changes*
        if state == self._prev_state:
            return

        if state == _DROGUE_P_DEPLOYING:
            self._drogue_p.set_deploying()

        elif state == _DROGUE_P_DEPLOYED:
            self._drogue_p.set_deployed(timestamp)

        elif state == _DROGUE_S_DEPLOYING:
            self._drogue_s.set_deploying()

        elif state == _DROGUE_S_DEPLOYED:
            self._drogue_s.set_deployed(timestamp)

        elif state == _MAIN_P_DEPLOYING:
            self._main_p.set_deploying()

        elif state == _MAIN_P_DEPLOYED:
            self._main_p.set_deployed(timestamp)

        elif state == _MAIN_S_DEPLOYING:
            self._main_s.set_deploying()

        elif state == _MAIN_S_DEPLOYED:
            self._main_s.set_deployed(timestamp)

        self._prev_state = state