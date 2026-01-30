import sys
from PyQt6.QtWidgets import QApplication
from Frontend.main_window import GroundStationWindow
import qdarktheme
from Backend.backend import list_serial_ports

def choose_port_from_console():
    ports = list_serial_ports()
    if not ports:
        print("No serial ports found.")
        return None

    print("Available serial ports:")
    for i, (dev, desc) in enumerate(ports):
        print(f"  [{i}] {dev} - {desc}")

    while True:
        sel = input("Select port index (or Enter to skip): ").strip()
        if sel == "":
            return None
        try:
            idx = int(sel)
            if 0 <= idx < len(ports):
                return ports[idx][0]
        except Exception:
            pass
        print("Invalid selection.")

def main():
    """Give the option to choose serial ports"""
    selected_port = choose_port_from_console()

    """Main entry point for the 0d Station application."""
    app = QApplication(sys.argv)
    
    # Apply modern dark theme
    app.setStyleSheet(qdarktheme.load_stylesheet()) #Not working for some reason
    
    window = GroundStationWindow()
    window.selected_port = selected_port
    window.show()
    
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
