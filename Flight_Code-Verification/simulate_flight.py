#!/usr/bin/env python3
"""
Offline flight-state simulator.

Reads a flight-data CSV, replays the critical state machine from the
Arduino firmware (excluding Kalman and airbrake logic), and produces:
  * A per-sample flight state trace.
  * An event timeline (launch detect, drogue/main firings, etc.).
  * A plot overlaying measured altitude with detected events.

Usage
-----
python simulate_flight.py /path/to/rocket.csv --output plot.png

Requirements
------------
pip install pandas matplotlib
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from enum import IntEnum, auto
from pathlib import Path
from typing import Deque, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd
from collections import deque
from itertools import cycle

# --- Constants copied from firmware (units preserved) ---
LAUNCH_THRESHOLD = 0.5
LANDED_THRESHOLD = -0.5
ALT_DIF_BUF_SIZE = 10

CHARGE_DELAY_MS = 500
BACKUP_DELAY_MS = 500

DROGUE_MAIN_MIN_ALT_M = 152
DROGUE_MAIN_MAX_ALT_M = 305

FALL_COUNTER_TRIGGER = 10


class FlightState(IntEnum):
    LAUNCH_PAD = auto()
    MOTOR_BURN = auto()
    GLIDING_ASCENT = auto()
    DROGUE_PRIMARY_DEPLOYING = auto()
    DROGUE_PRIMARY_DEPLOYED = auto()
    DROGUE_SECONDARY_DEPLOYING = auto()
    DROGUE_SECONDARY_DEPLOYED = auto()
    MAIN_PRIMARY_DEPLOYING = auto()
    MAIN_PRIMARY_DEPLOYED = auto()
    MAIN_SECONDARY_DEPLOYING = auto()
    MAIN_SECONDARY_DEPLOYED = auto()
    LANDED = auto()


STATE_LABELS = {state.value: state.name for state in FlightState}


@dataclass
class Event:
    time_s: float
    altitude_m: float
    kind: str
    state: FlightState
    note: str


def rolling_avg(buffer: Deque[float]) -> float:
    values = list(buffer)
    if len(values) < 3:
        return sum(values) / max(1, len(values))
    largest = max(values)
    smallest = min(values)
    return (sum(values) - largest - smallest) / (len(values) - 2)


def ensure_time_column(df: pd.DataFrame, time_col: Optional[str], sample_period_ms: Optional[float]) -> Tuple[pd.Series, str]:
    if time_col and time_col in df.columns:
        return df[time_col], time_col
    if time_col:
        print(f"[warn] Column '{time_col}' not found; synthesizing time column.")
    if sample_period_ms is None:
        raise ValueError("No time column found. Provide --time-column or --sample-period-ms.")
    synthesized = pd.Series(df.index * sample_period_ms, name="time_ms")
    df.insert(0, synthesized.name, synthesized)
    return synthesized, synthesized.name


def simulate_states(df: pd.DataFrame, alt_col: str, time_col: str) -> Tuple[List[FlightState], List[Event]]:
    if alt_col not in df.columns:
        raise ValueError(f"Altitude column '{alt_col}' not present in CSV.")

    start_alt = float(df[alt_col].iloc[0])
    pre_alt = start_alt
    alt_buf: Deque[float] = deque([0.0] * ALT_DIF_BUF_SIZE, maxlen=ALT_DIF_BUF_SIZE)

    state = FlightState.LAUNCH_PAD
    states: List[FlightState] = []
    events: List[Event] = []

    # timers
    launch_start_time = None
    drogue_primary_start = None
    drogue_primary_end = None
    drogue_secondary_start = None
    main_primary_start = None
    main_primary_end = None
    main_secondary_start = None

    fall_counter = 0
    fall_counter1 = 0

    for row in df.itertuples(index=False):
        time_ms = getattr(row, time_col)
        alt = getattr(row, alt_col)
        alt_buf.append(alt - pre_alt)
        avg_alt_dif = rolling_avg(alt_buf)
        time_s = time_ms / 1000.0

        def record(kind: str, note: str) -> None:
            events.append(Event(time_s=time_s, altitude_m=alt, kind=kind, state=state, note=note))

        # state machine logic
        if state == FlightState.LAUNCH_PAD:
            if avg_alt_dif > LAUNCH_THRESHOLD:
                launch_start_time = time_ms
                record("state", "Launch detected -> MOTOR_BURN")
                state = FlightState.MOTOR_BURN

        elif state == FlightState.MOTOR_BURN:
            if launch_start_time is not None and time_ms - launch_start_time > 5000:
                record("state", "Motor burn complete -> GLIDING_ASCENT")
                state = FlightState.GLIDING_ASCENT

        elif state == FlightState.GLIDING_ASCENT:
            if avg_alt_dif < 0:
                drogue_primary_start = time_ms
                record("pyro", "Fire drogue primary")
                state = FlightState.DROGUE_PRIMARY_DEPLOYING
            elif pre_alt - alt > 0.1:
                fall_counter += 1
            else:
                fall_counter = 0

        elif state == FlightState.DROGUE_PRIMARY_DEPLOYING:
            if drogue_primary_start is not None and time_ms - drogue_primary_start >= CHARGE_DELAY_MS:
                drogue_primary_end = time_ms
                record("state", "Drogue primary deployed")
                state = FlightState.DROGUE_PRIMARY_DEPLOYED

        elif state == FlightState.DROGUE_PRIMARY_DEPLOYED:
            if drogue_primary_end is not None and time_ms - drogue_primary_end >= BACKUP_DELAY_MS:
                drogue_secondary_start = time_ms
                record("pyro", "Fire drogue secondary")
                state = FlightState.DROGUE_SECONDARY_DEPLOYING

        elif state == FlightState.DROGUE_SECONDARY_DEPLOYING:
            if drogue_secondary_start is not None and time_ms - drogue_secondary_start >= CHARGE_DELAY_MS:
                record("state", "Drogue secondary deployed")
                state = FlightState.DROGUE_SECONDARY_DEPLOYED

        elif state == FlightState.DROGUE_SECONDARY_DEPLOYED:
            altitude_agl = alt - start_alt
            if (
                DROGUE_MAIN_MIN_ALT_M <= altitude_agl <= DROGUE_MAIN_MAX_ALT_M
                and pre_alt - alt > 1
                and fall_counter1 >= FALL_COUNTER_TRIGGER
            ):
                main_primary_start = time_ms
                record("pyro", "Fire main primary")
                state = FlightState.MAIN_PRIMARY_DEPLOYING
                fall_counter1 = 0
            elif pre_alt - alt > 0.1:
                fall_counter1 += 1
            else:
                fall_counter1 = 0

        elif state == FlightState.MAIN_PRIMARY_DEPLOYING:
            if main_primary_start is not None and time_ms - main_primary_start >= CHARGE_DELAY_MS:
                main_primary_end = time_ms
                record("state", "Main primary deployed")
                state = FlightState.MAIN_PRIMARY_DEPLOYED

        elif state == FlightState.MAIN_PRIMARY_DEPLOYED:
            if main_primary_end is not None and time_ms - main_primary_end >= BACKUP_DELAY_MS:
                main_secondary_start = time_ms
                record("pyro", "Fire main secondary")
                state = FlightState.MAIN_SECONDARY_DEPLOYING

        elif state == FlightState.MAIN_SECONDARY_DEPLOYING:
            if main_secondary_start is not None and time_ms - main_secondary_start >= CHARGE_DELAY_MS:
                record("state", "Main secondary deployed")
                state = FlightState.MAIN_SECONDARY_DEPLOYED

        elif state == FlightState.MAIN_SECONDARY_DEPLOYED:
            if avg_alt_dif > LANDED_THRESHOLD:
                record("state", "Landed")
                state = FlightState.LANDED

        # LANDED -> no further transitions

        states.append(state)
        pre_alt = alt

    return states, events


def plot_results(df: pd.DataFrame, time_col: str, alt_col: str, states: Iterable[FlightState], events: List[Event], output_path: Optional[Path]) -> None:
    time_s = df[time_col] / 1000.0
    altitude = df[alt_col]

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.plot(time_s, altitude, label="Measured Altitude", color="#1f77b4")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Altitude (m)")

    if events:
        marker_options = ["o", "s", "^", "D", "P", "X", "*", "v", ">", "<"]
        color_options = plt.cm.tab20.colors
        style_cycle = cycle([(m, c) for m in marker_options for c in color_options])
        note_styles: dict[str, tuple[str, str]] = {}
        for ev in events:
            if ev.note not in note_styles:
                note_styles[ev.note] = next(style_cycle)

        for note, (marker, color) in note_styles.items():
            xs = [ev.time_s for ev in events if ev.note == note]
            ys = [ev.altitude_m for ev in events if ev.note == note]
            ax.scatter(xs, ys, marker=marker, color=color, label=note, s=55, alpha=0.9)

    ax.legend(loc="best")

    fig.tight_layout()
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300)
        print(f"Plot saved to: {output_path}")
    else:
        plt.show()
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Replay rocket firmware state machine on recorded flight data.")
    parser.add_argument("csv_path", type=Path, help="CSV file containing at least altitude and time columns.")
    parser.add_argument("--alt-column", default="Alt", help="Column name for altitude (default: Alt).")
    parser.add_argument("--time-column", default="Time", help="Column name for timestamp in ms (default: Time).")
    parser.add_argument("--sample-period-ms", type=float, default=None, help="Sample period (ms) if no time column is present.")
    parser.add_argument("--output", type=Path, default=None, help="Path for the generated plot (PNG).")
    parser.add_argument("--events-csv", type=Path, default=None, help="Optional path to write the detected events CSV.")
    args = parser.parse_args()

    if not args.csv_path.exists():
        raise FileNotFoundError(args.csv_path)

    df = pd.read_csv(args.csv_path)
    time_series, time_col = ensure_time_column(df, args.time_column, args.sample_period_ms)

    states, events = simulate_states(df, args.alt_column, time_col)
    df["SimState"] = [state.name for state in states]

    plot_results(df, time_col, args.alt_column, states, events, args.output)

    if args.events_csv:
        events_df = pd.DataFrame(
            {
                "time_s": [ev.time_s for ev in events],
                "altitude_m": [ev.altitude_m for ev in events],
                "kind": [ev.kind for ev in events],
                "state": [ev.state.name for ev in events],
                "note": [ev.note for ev in events],
            }
        )
        events_df.to_csv(args.events_csv, index=False)
        print(f"Events saved to: {args.events_csv}")

    print("\nDetected events:")
    for ev in events:
        print(f"{ev.time_s:8.2f}s | {ev.altitude_m:8.2f} m | {ev.kind.upper():5s} | {ev.note}")


if __name__ == "__main__":
    main()

