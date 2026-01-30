# test_temp.py
import sys, time, random
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import QTimer
from Frontend.temp_graph import TempGraph

app = QApplication(sys.argv)
w = TempGraph()
w.show()

t0 = time.time()
temp = 25.0

def tick():
    global temp
    t = time.time() - t0
    temp += 0.02 * (random.random() - 0.5)  # small noise around baseline
    w.update_data(t, temp)

timer = QTimer()
timer.timeout.connect(tick)
timer.start(100)  # 10 Hz
sys.exit(app.exec())
