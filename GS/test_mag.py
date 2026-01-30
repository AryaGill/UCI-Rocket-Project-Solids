# test_mag.py
import sys, time, random
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import QTimer
from Frontend.mag_graph import MagGraph

app = QApplication(sys.argv)
w = MagGraph()
w.show()

t0 = time.time()

def tick():
    t = time.time() - t0
    # fake 3-axis sensor-ish values
    x = 50 * (random.random() - 0.5)
    y = 50 * (random.random() - 0.5)
    z = 50 * (random.random() - 0.5)
    w.update_data(t, x, y, z)

timer = QTimer()
timer.timeout.connect(tick)
timer.start(50)  # 20 Hz
sys.exit(app.exec())
