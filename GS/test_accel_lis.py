# test_accel_lis.py
import sys, time, random
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import QTimer
from Frontend.accel_graphlis import AccelGraphLIS

app = QApplication(sys.argv)
w = AccelGraphLIS()
w.show()

t0 = time.time()

def tick():
    t = time.time() - t0
    ax = 2.0 * (random.random() - 0.5)  # g-ish
    ay = 2.0 * (random.random() - 0.5)
    az = 1.0 + 0.2 * (random.random() - 0.5)  # bias around 1g
    w.update_data(t, ax, ay, az)

timer = QTimer()
timer.timeout.connect(tick)
timer.start(20)  # 50 Hz
sys.exit(app.exec())
