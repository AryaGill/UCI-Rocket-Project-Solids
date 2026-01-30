from PyQt6.QtWidgets import QWidget, QFrame, QLabel, QHBoxLayout
from PyQt6.QtCore import Qt

class StatusIndicator(QWidget):
    """A circular status indicator with text: green = in-flight, red = inactive."""

    def __init__(self, diameter=30, parent=None):
        super().__init__(parent)

        self.diameter = diameter

        # Layout for circle + label
        layout = QHBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(10)  # space between circle and text
        layout.setAlignment(Qt.AlignmentFlag.AlignRight)  # stick to right
        self.setLayout(layout)

        # Circle indicator
        self.circle = QFrame()
        self.circle.setFixedSize(self.diameter, self.diameter)
        self.circle.setFrameShape(QFrame.Shape.NoFrame)
        self.circle.setStyleSheet(f"border-radius: {self.diameter // 2}px; background-color: red;")
        layout.addWidget(self.circle)

        # Status text
        self.label = QLabel("Inactive")
        layout.addWidget(self.label)

        # Default state
        self.set_status(False)

    def set_status(self, in_flight: bool):
        """Update the indicator color and text."""
        color = "green" if in_flight else "red"
        text = "Active" if in_flight else "Inactive"
        self.circle.setStyleSheet(f"border-radius: {self.diameter // 2}px; background-color: {color};")
        self.label.setText(text)
