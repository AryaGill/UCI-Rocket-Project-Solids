import sys
from PyQt6.QtWidgets import QApplication
# import qdarktheme
from Backend.backend import select_serial_port
from Frontend.main_window import GroundStationWindow


def main():
    """Main entry point for the Ground Station application."""
    # Get port selection before starting GUI
    selected_port = select_serial_port()

    app = QApplication(sys.argv)

    # Apply modern dark theme
    # app.setStyleSheet(qdarktheme.load_stylesheet())

    # Pass port directly to window constructor
    window = GroundStationWindow(port=selected_port)
    window.show()

    sys.exit(app.exec())


if __name__ == '__main__':
    main()