from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QGridLayout,
                             QPushButton, QLabel, QMessageBox)
from PyQt6.QtCore import Qt, pyqtSignal

class CameraPanel(QDialog):
    """
    Popup dialog for pyrotechnic charge control.
    Sends RF commands over serial to fire parachute deployment charges.
    """
    
    command_signal = pyqtSignal(str)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Camera Control")
        self.setModal(False)  # Allow interaction with cam2 window
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
        
        # cam1s (Row 0)
        cam1_label = QLabel("Camera 1")
        cam1_label.setStyleSheet("""
            QLabel {
                color: #00d4ff;
                font-weight: bold;
                font-size: 13px;
                padding: 5px;
            }
        """)
        button_layout.addWidget(cam1_label, 0, 0, 1, 2, Qt.AlignmentFlag.AlignCenter)
        
        self.cam1_p_btn = self.create_pyro_button("Camera 1 ON", "#007bff")
        self.cam1_p_btn.clicked.connect(lambda: self.send_command("CAM1ON", 1))
        button_layout.addWidget(self.cam1_p_btn, 1, 0)
        
        self.cam1_s_btn = self.create_pyro_button("Camera 1 OFF", "#007bff")
        self.cam1_s_btn.clicked.connect(lambda: self.send_command("CAM1OFF", 1))
        button_layout.addWidget(self.cam1_s_btn, 1, 1)
        
        # cam2 (Row 2)
        cam2_label = QLabel("Camera 2")
        cam2_label.setStyleSheet("""
            QLabel {
                color: #28dc5e;
                font-weight: bold;
                font-size: 13px;
                padding: 5px;
            }
        """)
        button_layout.addWidget(cam2_label, 2, 0, 1, 2, Qt.AlignmentFlag.AlignCenter)
        
        self.cam2_p_btn = self.create_pyro_button("Camera 2 ON", "#9c60f6")
        self.cam2_p_btn.clicked.connect(lambda: self.send_command("CAM2ON", 2))
        button_layout.addWidget(self.cam2_p_btn, 3, 0)
        
        self.cam2_s_btn = self.create_pyro_button("Camera 2 OFF", "#9c60f6")
        self.cam2_s_btn.clicked.connect(lambda: self.send_command("CAM2OFF", 2))
        button_layout.addWidget(self.cam2_s_btn, 3, 1)

        #Dual cams (Row 3)
        dual_label = QLabel("Both Cameras")
        dual_label.setStyleSheet("""
            QLabel {
                color: #28dc5e;
                font-weight: bold;
                font-size: 13px;
                padding: 5px;
            }
        """)
        button_layout.addWidget(dual_label, 4, 0, 1, 2, Qt.AlignmentFlag.AlignCenter)
        
        self.dual_p_btn = self.create_pyro_button("Both Cameras ON", "#ff36e4")
        self.dual_p_btn.clicked.connect(lambda: self.send_command("CAM1ON", 1))
        self.dual_p_btn.clicked.connect(lambda: self.send_command("CAM2ON", 2))
        button_layout.addWidget(self.dual_p_btn, 5, 0)
        
        self.dual_s_btn = self.create_pyro_button("Both Cameras OFF", "#ff36e4")
        self.dual_s_btn.clicked.connect(lambda: self.send_command("CAM1OFF", 1))
        self.dual_s_btn.clicked.connect(lambda: self.send_command("CAM2OFF", 2))
        button_layout.addWidget(self.dual_s_btn, 5, 1)
        
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
        self.command_signal.emit(command)
