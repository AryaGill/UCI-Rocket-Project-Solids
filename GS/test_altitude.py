# test_altitude.py
import sys, time, random
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import QTimer
from Frontend.altitude_graph import AltitudeGraph

app = QApplication(sys.argv)
w = AltitudeGraph()
w.show()

t0 = time.time()
alt = 0.0

def tick():
    global alt
    t = time.time() - t0
    alt += 0.8 + 0.2 * (random.random() - 0.5)  # slow climb
    w.update_data(t, alt)

timer = QTimer()
timer.timeout.connect(tick)
timer.start(50)  # 20 Hz
sys.exit(app.exec())
