from PyQt6.QtWidgets import QFrame, QGridLayout, QPushButton, QSizePolicy
from PyQt6.QtCore import pyqtSignal

class PyroDroguePanel(QFrame):
    arm_clicked = pyqtSignal()
    fire_clicked = pyqtSignal()
    disarm_clicked = pyqtSignal()
    abort_clicked = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)

        self.setVisible(False)
        self.setStyleSheet("""
            QFrame {
                background-color: #1e1e1e;
                border-top: 2px solid #404040;
            }
        """)

        layout = QGridLayout(self)
        layout.setContentsMargins(20, 10, 20, 10)
        layout.setHorizontalSpacing(20)
        layout.setVerticalSpacing(12)

        self.arm_btn = QPushButton("Main Primary")
        self.fire_btn = QPushButton("Main Secondary")
        self.disarm_btn = QPushButton("Drogue Primary")
        self.abort_btn = QPushButton("Drogue Secondary")

        layout.addWidget(self.arm_btn,    0, 0)
        layout.addWidget(self.fire_btn,   0, 1)
        layout.addWidget(self.disarm_btn, 1, 0)
        layout.addWidget(self.abort_btn,  1, 1)

        for btn in (
            self.arm_btn,
            self.fire_btn,
            self.disarm_btn,
            self.abort_btn
        ):
            btn.setMinimumHeight(36)
            btn.setSizePolicy(QSizePolicy.Policy.Expanding,
                              QSizePolicy.Policy.Fixed)

        layout.addWidget(self.arm_btn)
        layout.addWidget(self.fire_btn)
        layout.addWidget(self.disarm_btn)
        layout.addWidget(self.abort_btn)

        # Signal wiring
        self.arm_btn.clicked.connect(self.arm_clicked)
        self.fire_btn.clicked.connect(self.fire_clicked)
        self.disarm_btn.clicked.connect(self.disarm_clicked)
        self.abort_btn.clicked.connect(self.abort_clicked)
