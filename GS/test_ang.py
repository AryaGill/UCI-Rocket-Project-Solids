# test_ang.py
import sys, time, random
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import QTimer
from Frontend.ang_graph import AngGraph

app = QApplication(sys.argv)
w = AngGraph()
w.show()

t0 = time.time()

def tick():
    t = time.time() - t0
    gx = 200 * (random.random() - 0.5)  # deg/s-ish
    gy = 200 * (random.random() - 0.5)
    gz = 200 * (random.random() - 0.5)
    w.update_data(t, gx, gy, gz)

timer = QTimer()
timer.timeout.connect(tick)
timer.start(50)  # 20 Hz
sys.exit(app.exec())
