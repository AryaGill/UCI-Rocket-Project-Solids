import sys
from PyQt6.QtWidgets import QApplication
from Frontend.main_window import GroundStationWindow
import qdarktheme


def main():
    """Main entry point for the Ground Station application."""
    app = QApplication(sys.argv)
    
    # Apply modern dark theme
    qdarktheme.setup_theme("dark")
    
    window = GroundStationWindow()
    window.show()
    
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
