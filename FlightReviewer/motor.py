import numpy as np
import pandas as pd

# ── Detection thresholds ───────────────────────────────────────────────────────
LAUNCH_THRESH_G   = 1.5   # g  — sustained accel above this = motor burning
BURNOUT_THRESH_G  = 0.5   # g  — sustained accel below this (after ignition) = burnout
CONFIRM_SAMPLES   = 2     # consecutive samples needed to confirm state change
SMOOTHING_WINDOW  = 5     # rolling-mean window (samples) applied before threshold check

MOTOR_CLASSES = [
    (0.3125,    "1/8A"), (0.625,   "1/4A"), (1.25,    "1/2A"),
    (2.5,       "A"),    (5.0,     "B"),     (10.0,    "C"),
    (20.0,      "D"),    (40.0,    "E"),     (80.0,    "F"),
    (160.0,     "G"),    (320.0,   "H"),     (640.0,   "I"),
    (1280.0,    "J"),    (2560.0,  "K"),     (5120.0,  "L"),
    (10240.0,   "M"),    (20480.0, "N"),     (40960.0, "O"),
]


def classify_motor(total_impulse_ns: float) -> str:
    for max_ns, letter in MOTOR_CLASSES:
        if total_impulse_ns <= max_ns:
            return letter
    return "P+"


def _smooth_g(accel_raw: np.ndarray) -> np.ndarray:
    """Smooth raw accel (m/s²) and convert to g."""
    smoothed = pd.Series(accel_raw).rolling(SMOOTHING_WINDOW, center=True).mean().to_numpy()
    return smoothed / 9.81


def detect_burn_start(time_ms: np.ndarray, z_accel_raw: np.ndarray):
    """Return (index, time_ms) of confirmed burn start, or (None, None)."""
    z_g = _smooth_g(z_accel_raw)
    consecutive = 0
    for i in range(len(z_g)):
        if z_g[i] > LAUNCH_THRESH_G:
            consecutive += 1
            if consecutive >= CONFIRM_SAMPLES:
                idx = i - CONFIRM_SAMPLES + 1
                return idx, time_ms[idx]
        else:
            consecutive = 0
    return None, None


def detect_burn_end(time_ms: np.ndarray, z_accel_raw: np.ndarray, burn_start_idx):
    """Return (index, time_ms) of confirmed burnout, or (None, None)."""
    if burn_start_idx is None:
        return None, None
    z_g = _smooth_g(z_accel_raw)
    consecutive = 0
    for i in range(burn_start_idx, len(z_g)):
        if z_g[i] < BURNOUT_THRESH_G:
            consecutive += 1
            if consecutive >= CONFIRM_SAMPLES:
                idx = i - CONFIRM_SAMPLES + 1
                return idx, time_ms[idx]
        else:
            consecutive = 0
    return None, None


def derive_thrust(accel_world_z: np.ndarray, mass_kg: float) -> np.ndarray:
    """F = m * a_world_z. Clamp negatives to 0 — motor only pushes."""
    return np.maximum(accel_world_z * mass_kg, 0.0)


def compute_motor_stats(time_ms: np.ndarray, thrust: np.ndarray,
                        accel_world_z: np.ndarray = None) -> dict:
    """
    Detect burn window using g-threshold logic (if accel_world_z provided),
    then integrate thrust over that window.
    time_ms is in milliseconds; all outputs use seconds.
    """
    if accel_world_z is not None:
        i0, _ = detect_burn_start(time_ms, accel_world_z)
        if i0 is None:
            return _empty_stats()
        i1, _ = detect_burn_end(time_ms, accel_world_z, i0)
        if i1 is None:
            return _empty_stats()
    else:
        nonzero = np.where(thrust > 0)[0]
        if len(nonzero) == 0:
            return _empty_stats()
        i0, i1 = int(nonzero[0]), int(nonzero[-1])

    time_s        = time_ms / 1000.0
    t_burn        = time_s[i0:i1 + 1]
    f_burn        = thrust[i0:i1 + 1]
    burn_time     = float(t_burn[-1] - t_burn[0])
    integrate      = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    total_impulse = float(integrate(f_burn, t_burn))
    avg_thrust    = total_impulse / burn_time if burn_time > 0 else 0.0
    peak_thrust   = float(np.max(f_burn))
    designation   = classify_motor(total_impulse)

    return {
        "total_impulse": total_impulse,
        "burn_time":     burn_time,
        "avg_thrust":    avg_thrust,
        "peak_thrust":   peak_thrust,
        "designation":   designation,
    }


def _empty_stats() -> dict:
    return {k: None for k in ("total_impulse", "burn_time", "avg_thrust", "peak_thrust", "designation")}
