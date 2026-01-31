from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QGridLayout,
                             QPushButton, QLabel, QMessageBox)
from PyQt6.QtCore import Qt, pyqtSignal

class PyroPanel(QDialog):
    """
    Popup dialog for pyrotechnic charge control.
    Sends RF commands over serial to fire parachute deployment charges.
    """
    
    command_signal = pyqtSignal(str)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Pyro Charges")
        self.setModal(False)  # Allow interaction with main window
        self.setMinimumWidth(400)
        self.setup_ui()
    
    def setup_ui(self):
        """Setup the pyro control interface."""
        layout = QVBoxLayout()
        
        # Warning label
        warning = QLabel("⚠️ WARNING: PYROTECHNIC DEVICES")
        warning.setStyleSheet("""
            QLabel {
                background-color: #ff6b35;
                color: #1e1e1e;
                font-weight: bold;
                font-size: 14px;
                padding: 10px;
                border-radius: 5px;
            }
        """)
        warning.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(warning)
        
        # Instructions
        info = QLabel("Click buttons to send fire commands via RF")
        info.setStyleSheet("color: #b0b0b0; padding: 5px;")
        info.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(info)
        
        # Button grid
        button_layout = QGridLayout()
        
        # Drogue charges (Row 0)
        drogue_label = QLabel("Drogue")
        drogue_label.setStyleSheet("""
            QLabel {
                color: #00d4ff;
                font-weight: bold;
                font-size: 13px;
                padding: 5px;
            }
        """)
        button_layout.addWidget(drogue_label, 0, 0, 1, 2, Qt.AlignmentFlag.AlignCenter)
        
        self.drogue_p_btn = self.create_pyro_button("Drogue Primary", "#ff6b35")
        self.drogue_p_btn.clicked.connect(lambda: self.send_command("Fire Drogue P"))
        button_layout.addWidget(self.drogue_p_btn, 1, 0)
        
        self.drogue_s_btn = self.create_pyro_button("Drogue Secondary", "#ff8555")
        self.drogue_s_btn.clicked.connect(lambda: self.send_command("Fire Drogue S"))
        button_layout.addWidget(self.drogue_s_btn, 1, 1)
        
        # Main charges (Row 2)
        main_label = QLabel("Main")
        main_label.setStyleSheet("""
            QLabel {
                color: #28dc5e;
                font-weight: bold;
                font-size: 13px;
                padding: 5px;
            }
        """)
        button_layout.addWidget(main_label, 2, 0, 1, 2, Qt.AlignmentFlag.AlignCenter)
        
        self.main_p_btn = self.create_pyro_button("Main Primary", "#28dc5e")
        self.main_p_btn.clicked.connect(lambda: self.send_command("Fire Main P"))
        button_layout.addWidget(self.main_p_btn, 3, 0)
        
        self.main_s_btn = self.create_pyro_button("Main Secondary", "#52e67d")
        self.main_s_btn.clicked.connect(lambda: self.send_command("Fire Main S"))
        button_layout.addWidget(self.main_s_btn, 3, 1)
        
        layout.addLayout(button_layout)
        
        # Close button
        close_btn = QPushButton("Close")
        close_btn.setStyleSheet("""
            QPushButton {
                background-color: #404040;
                color: #b0b0b0;
                border: none;
                padding: 8px;
                font-weight: bold;
                border-radius: 3px;
            }
            QPushButton:hover {
                background-color: #505050;
            }
        """)
        close_btn.clicked.connect(self.close)
        layout.addWidget(close_btn)
        
        self.setLayout(layout)
    
    def create_pyro_button(self, text, color):
        """Create a styled pyro button."""
        btn = QPushButton(text)
        btn.setMinimumHeight(50)
        btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {color};
                color: #1e1e1e;
                border: 2px solid {color};
                padding: 10px;
                font-weight: bold;
                font-size: 12px;
                border-radius: 5px;
            }}
            QPushButton:hover {{
                background-color: #1e1e1e;
                color: {color};
            }}
            QPushButton:pressed {{
                background-color: #0d0d0d;
            }}
        """)
        return btn
    
    def send_command(self, command):
        """Confirm and send pyro command."""
        # Confirmation dialog
        reply = QMessageBox.question(
            self,
            "Confirm Pyro Command",
            f"Send command: '{command}'?\n\nThis will fire a pyrotechnic charge!",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No
        )
        
        if reply == QMessageBox.StandardButton.Yes:
            self.command_signal.emit(command)
