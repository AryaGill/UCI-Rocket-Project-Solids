import time, random
from PyQt6.QtCore import QThread, pyqtSignal

class FakeStreamer(QThread):
    new_data = pyqtSignal(dict)
    finished = pyqtSignal()

    def __init__(self, hz=20, parent=None):
        super().__init__(parent)
        self.hz = hz
        self.is_running = True

    def run(self):
        t = 0.0
        dt = 1.0 / self.hz
        while self.is_running:
            t += dt
            self.new_data.emit({
                "Time": t,
                "Alt": 100 + 5 * (random.random() - 0.5),
                "Temp": 25 + 2 * (random.random() - 0.5),
                "MagX": 50 * (random.random() - 0.5),
                "MagY": 50 * (random.random() - 0.5),
                "MagZ": 50 * (random.random() - 0.5),
            })
            time.sleep(dt)
        self.finished.emit()

    def stop(self):
        self.is_running = False
