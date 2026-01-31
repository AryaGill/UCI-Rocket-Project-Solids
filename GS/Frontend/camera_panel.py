from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout,
                             QPushButton, QLabel)
from PyQt6.QtCore import Qt, pyqtSignal

class CameraPanel(QDialog):
    """
    Popup dialog for camera control.
    Sends RF commands over serial to control camera.
    """
    
    command_signal = pyqtSignal(str)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Camera Control")
        self.setModal(False)  # Allow interaction with main window
        self.setMinimumWidth(350)
        self.setup_ui()
    
    def setup_ui(self):
        """Setup the camera control interface."""
        layout = QVBoxLayout()
        
        # Header label
        header = QLabel("📹 Camera Control")
        header.setStyleSheet("""
            QLabel {
                background-color: #fffb00;
                color: #1e1e1e;
                font-weight: bold;
                font-size: 14px;
                padding: 10px;
                border-radius: 5px;
            }
        """)
        header.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(header)
        
        # Instructions
        info = QLabel("Control camera via RF commands")
        info.setStyleSheet("color: #b0b0b0; padding: 5px;")
        info.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(info)
        
        # Camera buttons
        camera_layout = QHBoxLayout()
        
        self.camera_on_btn = QPushButton("Camera ON")
        self.camera_on_btn.setMinimumHeight(50)
        self.camera_on_btn.setStyleSheet("""
            QPushButton {
                background-color: #28dc5e;
                color: #1e1e1e;
                border: 2px solid #28dc5e;
                padding: 10px;
                font-weight: bold;
                font-size: 13px;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #1e1e1e;
                color: #28dc5e;
            }
            QPushButton:pressed {
                background-color: #0d0d0d;
            }
        """)
        self.camera_on_btn.clicked.connect(lambda: self.send_command("ON"))
        camera_layout.addWidget(self.camera_on_btn)
        
        self.camera_off_btn = QPushButton("Camera OFF")
        self.camera_off_btn.setMinimumHeight(50)
        self.camera_off_btn.setStyleSheet("""
            QPushButton {
                background-color: #ff6b35;
                color: #1e1e1e;
                border: 2px solid #ff6b35;
                padding: 10px;
                font-weight: bold;
                font-size: 13px;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #1e1e1e;
                color: #ff6b35;
            }
            QPushButton:pressed {
                background-color: #0d0d0d;
            }
        """)
        self.camera_off_btn.clicked.connect(lambda: self.send_command("OFF"))
        camera_layout.addWidget(self.camera_off_btn)
        
        layout.addLayout(camera_layout)
        
        # Status label
        self.status_label = QLabel("")
        self.status_label.setStyleSheet("""
            QLabel {
                color: #00d4ff;
                padding: 5px;
                font-size: 11px;
            }
        """)
        self.status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.status_label)
        
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
                margin-top: 10px;
            }
            QPushButton:hover {
                background-color: #505050;
            }
        """)
        close_btn.clicked.connect(self.close)
        layout.addWidget(close_btn)
        
        self.setLayout(layout)
    
    def send_command(self, command):
        """Send camera command."""
        self.command_signal.emit(command)
        self.status_label.setText(f"Sent: {command}")
