from PyQt6.QtWidgets import QWidget, QVBoxLayout, QLabel
from PyQt6.QtCore import Qt

class FlightStateDisplay(QWidget):
    """
    Widget to display current flight state from telemetry.
    """
    
    # Flight state mapping
    FLIGHT_STATES = {
        0: "DISARMED",
        1: "LAUNCH PAD",
        2: "MOTOR BURN",
        3: "GLIDING ASCENT",
        4: "DROGUE P DEPLOYING",
        5: "DROGUE P DEPLOYED",
        6: "DROGUE S DEPLOYING",
        7: "DROGUE S DEPLOYED",
        8: "MAIN P DEPLOYING",
        9: "MAIN P DEPLOYED",
        10: "MAIN S DEPLOYING",
        11: "MAIN S DEPLOYED",
        12: "LANDED"
    }
    
    # Color mapping for different states
    STATE_COLORS = {
        0: "#ff6b35",  # Orange - Disarmed
        1: "#fffb00",  # Yellow - Launch Pad
        2: "#ff6b35",  # Orange - Motor Burn
        3: "#00d4ff",  # Cyan - Gliding
        4: "#ff6b35",  # Orange - Drogue P Deploying
        5: "#28dc5e",  # Green - Drogue P Deployed
        6: "#ff6b35",  # Orange - Drogue S Deploying
        7: "#28dc5e",  # Green - Drogue S Deployed
        8: "#ff6b35",  # Orange - Main P Deploying
        9: "#28dc5e",  # Green - Main P Deployed
        10: "#ff6b35", # Orange - Main S Deploying
        11: "#28dc5e", # Green - Main S Deployed
        12: "#00d4ff"  # Cyan - Landed
    }
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.current_state = None
        self.setup_ui()
    
    def setup_ui(self):
        """Setup the flight state display interface."""
        layout = QVBoxLayout()
        
        # Title label
        title = QLabel("Flight State")
        title.setStyleSheet("""
            QLabel {
                color: #b0b0b0;
                font-weight: bold;
                font-size: 12px;
                padding: 5px;
            }
        """)
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)
        
        # State display label
        self.state_label = QLabel("UNKNOWN")
        self.state_label.setStyleSheet("""
            QLabel {
                background-color: #404040;
                color: #b0b0b0;
                font-weight: bold;
                font-size: 16px;
                padding: 15px;
                border-radius: 5px;
                border: 2px solid #404040;
            }
        """)
        self.state_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.state_label.setMinimumHeight(60)
        layout.addWidget(self.state_label)
        
        layout.setContentsMargins(5, 5, 5, 5)
        self.setLayout(layout)
    
    def update_state(self, state_value):
        """
        Update the displayed flight state.
        
        Args:
            state_value: Integer or string representing flight state (0-12)
        """
        try:
            # Convert to int if string
            if isinstance(state_value, str):
                state_value = int(state_value)
            
            # Get state name
            state_name = self.FLIGHT_STATES.get(state_value, "UNKNOWN")
            
            # Get state color
            state_color = self.STATE_COLORS.get(state_value, "#404040")
            
            # Update label
            self.state_label.setText(state_name)
            self.state_label.setStyleSheet(f"""
                QLabel {{
                    background-color: {state_color};
                    color: #1e1e1e;
                    font-weight: bold;
                    font-size: 16px;
                    padding: 15px;
                    border-radius: 5px;
                    border: 2px solid {state_color};
                }}
            """)
            
            self.current_state = state_value
            
        except (ValueError, TypeError):
            self.state_label.setText("INVALID")
            self.state_label.setStyleSheet("""
                QLabel {
                    background-color: #ff0000;
                    color: #1e1e1e;
                    font-weight: bold;
                    font-size: 16px;
                    padding: 15px;
                    border-radius: 5px;
                    border: 2px solid #ff0000;
                }
            """)
