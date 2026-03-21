from __future__ import annotations

from datetime import datetime

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QTextEdit, QPushButton, QLineEdit,
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QTextCursor, QColor, QTextCharFormat, QFont


class CommandEchoPanel(QWidget):
    """
    Scrollable command-echo console that logs every command sent *to* the
    flight computer and every status/acknowledgement message received *from*
    the backend.

    Public API
    ----------
    log_sent(cmd)      – call whenever a command is written to serial
    log_received(msg)  – call to append a status / telemetry note
    log_info(msg)      – neutral info lines (grey)
    clear()            – wipe the log

    The panel also exposes a manual-entry bar so operators can type
    arbitrary commands directly (emits ``command_entered`` signal).
    """

    command_entered: pyqtSignal = pyqtSignal(str)

    # Colour palette
    _COL_SENT     = "#fffb00"   # yellow  – commands we sent
    _COL_RECEIVED = "#28dc5e"   # green   – ack / telemetry notes
    _COL_INFO     = "#808080"   # grey    – neutral info
    _COL_ERROR    = "#ff4040"   # red     – errors / warnings
    _COL_TIMESTAMP= "#505050"   # dim grey for timestamps

    def __init__(self, parent=None):
        super().__init__(parent)
        self._build_ui()

    # ── UI construction ──────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(6, 2, 6, 2)
        root.setSpacing(2)

        # ── header row ───────────────────────────────────────────────────────
        header_row = QHBoxLayout()
        header_row.setSpacing(1)

        title = QLabel("FC Echo  (received)")
        title.setStyleSheet(
            "color: #b0b0b0; font-size: 10px; font-weight: bold;"
        )
        header_row.addWidget(title)
        header_row.addStretch()

        clear_btn = QPushButton("Clear")
        clear_btn.setFixedHeight(18)
        clear_btn.setFixedWidth(46)
        clear_btn.setStyleSheet("""
            QPushButton {
                background-color: #333333;
                color: #808080;
                border: 1px solid #404040;
                border-radius: 3px;
                font-size: 9px;
                padding: 0 4px;
            }
            QPushButton:hover { background-color: #444444; color: #b0b0b0; }
            QPushButton:pressed { background-color: #222222; }
        """)
        clear_btn.clicked.connect(self.clear)
        header_row.addWidget(clear_btn)

        root.addLayout(header_row)

        # ── log area ─────────────────────────────────────────────────────────
        self._log = QTextEdit()
        self._log.setReadOnly(True)
        self._log.setMinimumHeight(100)
        self._log.setMaximumHeight(100)
        self._log.setFont(QFont("Courier New", 9))
        self._log.setStyleSheet("""
            QTextEdit {
                background-color: #181818;
                color: #c0c0c0;
                border: 1px solid #383838;
                border-radius: 4px;
                padding: 4px;
            }
            QScrollBar:vertical {
                background: #1e1e1e; width: 8px; margin: 0;
            }
            QScrollBar::handle:vertical {
                background: #404040; border-radius: 4px; min-height: 20px;
            }
        """)
        root.addWidget(self._log)

        self.setFixedHeight(130)

    # ── public API ───────────────────────────────────────────────────────────

    def log_received(self, msg: str):
        """Log an ack / status message from the FC (green)."""
        self._append(f"◀ {msg}", self._COL_RECEIVED)

    def log_info(self, msg: str):
        """Log a neutral info line (grey)."""
        self._append(f"  {msg}", self._COL_INFO)

    def log_error(self, msg: str):
        """Log a warning / error line (red)."""
        self._append(f"✖ {msg}", self._COL_ERROR)

    def clear(self):
        self._log.clear()

    # ── internals ────────────────────────────────────────────────────────────

    def _timestamp(self) -> str:
        return datetime.now().strftime("%H:%M:%S")

    def _append(self, text: str, hex_color: str):
        cursor = self._log.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)

        # Timestamp in dim grey
        ts_fmt = QTextCharFormat()
        ts_fmt.setForeground(QColor(self._COL_TIMESTAMP))
        cursor.insertText(f"[{self._timestamp()}] ", ts_fmt)

        # Message in accent colour
        msg_fmt = QTextCharFormat()
        msg_fmt.setForeground(QColor(hex_color))
        cursor.insertText(text + "\n", msg_fmt)

        # Auto-scroll to bottom
        self._log.setTextCursor(cursor)
        self._log.ensureCursorVisible()

    def _on_send(self):
        cmd = self._entry.text().strip()
        if not cmd:
            return
        self._entry.clear()
        self.command_entered.emit(cmd)   # let the parent handle sending