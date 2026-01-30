# test_accel_lsm.py
import sys, time, random
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import QTimer
from Frontend.accel_graphlsm import AccelGraphLSM

app = QApplication(sys.argv)
w = AccelGraphLSM()
w.show()

t0 = time.time()

def tick():
    t = time.time() - t0
    ax = 16.0 * (random.random() - 0.5)  # larger range
    ay = 16.0 * (random.random() - 0.5)
    az = 9.8 + 0.5 * (random.random() - 0.5)
    w.update_data(t, ax, ay, az)

timer = QTimer()
timer.timeout.connect(tick)
timer.start(20)  # 50 Hz
sys.exit(app.exec())
