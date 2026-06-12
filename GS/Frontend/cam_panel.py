from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QGridLayout,
                             QPushButton, QLabel, QMessageBox)
from PyQt6.QtCore import Qt, pyqtSignal

class CamPanel(QDialog):
    """
    Popup dialog for pyrotechnic charge control.
    Sends RF commands over serial to fire parachute deployment charges.
    """
    
    command_signal = pyqtSignal(str)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Camera Control")
        self.setModal(False)  # Allow interaction with main window
        self.setMinimumWidth(400)
        self.setup_ui()
    
    def setup_ui(self):
        """Setup the cam control interface."""
        layout = QVBoxLayout()
        
        
        # Instructions
        info = QLabel("Click buttons to send fire commands via RF")
        info.setStyleSheet("color: #b0b0b0; padding: 5px;")
        info.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(info)
        
        # Button grid
        button_layout = QGridLayout()
        
        # Drogue charges (Row 0)
        drogue_label = QLabel("Camera 1")
        drogue_label.setStyleSheet("""
            QLabel {
                color: #00d4ff;
                font-weight: bold;
                font-size: 13px;
                padding: 5px;
            }
        """)
        button_layout.addWidget(drogue_label, 0, 0, 1, 2, Qt.AlignmentFlag.AlignCenter)
        
        self.drogue_p_btn = self.create_pyro_button("Camera 1 ON", "#007bff")
        self.drogue_p_btn.clicked.connect(lambda: self.send_command("CAM1ON"), 1)
        button_layout.addWidget(self.drogue_p_btn, 1, 0)
        
        self.drogue_s_btn = self.create_pyro_button("Camera 1 OFF", "#007bff")
        self.drogue_s_btn.clicked.connect(lambda: self.send_command("CAM1OFF"), 1)
        button_layout.addWidget(self.drogue_s_btn, 1, 1)
        
        # Main charges (Row 2)
        main_label = QLabel("Camera 2")
        main_label.setStyleSheet("""
            QLabel {
                color: #28dc5e;
                font-weight: bold;
                font-size: 13px;
                padding: 5px;
            }
        """)
        button_layout.addWidget(main_label, 2, 0, 1, 2, Qt.AlignmentFlag.AlignCenter)
        
        self.main_p_btn = self.create_pyro_button("Camera 2 ON", "#ff36e4")
        self.main_p_btn.clicked.connect(lambda: self.send_command("CAM2ON"), 2)
        button_layout.addWidget(self.main_p_btn, 3, 0)
        
        self.main_s_btn = self.create_pyro_button("Camera 2 OFF", "#ff36e4")
        self.main_s_btn.clicked.connect(lambda: self.send_command("CAM2OFF"), 2)
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
        """Create a styled cam button."""
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
    
    def send_command(self, command, num):
        """Confirm and send camera command."""
        # Confirmation dialog
        reply = QMessageBox.question(
            self,
            "Confirm Pyro Command",
            f"Send command: '{command}'?\n\nThis will affect camera {num}",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No
        )
        
        if reply == QMessageBox.StandardButton.Yes:
            self.command_signal.emit(command)
