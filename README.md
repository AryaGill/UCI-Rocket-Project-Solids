 # Flight Code Verification Toolkit

Tools for replaying and visualizing high-power rocket flights using recorded telemetry logs. This branch focuses on validating stage detection logic, cleaning raw datasets, and producing plots for pre-flight reviews.

## Requirements

- Python 3.10+
- `pip install pandas matplotlib`

## Scripts

| Script | Purpose | Key Flags |
| --- | --- | --- |
| `Flight_Code-Verification/simulate_flight.py` | Replays the Arduino flight state machine (sans Kalman + airbrakes) on a CSV, detects launch → recovery events, and plots altitude with event markers. | `--alt-column`, `--time-column`, `--output`, `--events-csv` |
| `Flight_Code-Verification/basic_plot.py` | Quick-look plotter for arbitrary CSV columns; handles duplicate headers. | `--x-axis`, `--output` |
| `Flight_Code-Verification/clean_flight_segment.py` | Extracts the flight containing the max altitude, trims to actual motion, and drops pre/post noise. | `--diff-threshold`, `--padding`, `--stable-window` |

## Typical Workflow

1. **Clean raw log** (optional):  
   ```bash
   python Flight_Code-Verification/clean_flight_segment.py raw.csv clean.csv \
     --diff-threshold 0.2 --padding 50 --stable-window 300
   ```

2. **Simulate flight computer logic**:  
   ```bash
   python Flight_Code-Verification/simulate_flight.py clean.csv \
     --alt-column Alt --time-column Time --events-csv events.csv
   ```
   - Shows the altitude plot with distinct markers for each event (launch detection, pyro fires, etc.).
   - Adds a `SimState` column to the CSV so you can diff against onboard state logging.

3. **Ad-hoc plotting** (e.g., compare accelerometer axes):  
   ```bash
   python Flight_Code-Verification/basic_plot.py clean.csv Accel_x Accel_y Accel_z --x-axis Time
   ```

## Notes

- All scripts default to displaying plots interactively; pass `--output path.png` to save.
- `simulate_flight.py` synthesizes a time column if the CSV lacks one (use `--sample-period-ms`).
- The cleaning script assumes duplicate header rows separate segments; only the peak-altitude segment is retained.
