#!/usr/bin/env python3
"""
Extract the primary flight segment from a CSV that contains multiple
concatenated data blocks separated by duplicate headers.

Strategy:
    1. Find the global peak altitude.
    2. Identify the header immediately preceding that peak.
    3. Keep rows from that header up to (but not including) the next header.
    4. Trim that block so it only includes the portion where altitude is
       actually changing (plus configurable padding before/after motion).

Usage:
    python clean_flight_segment.py input.csv output.csv
        [--diff-threshold 0.2] [--padding 200]

This keeps the flight containing the highest altitude, removes other
sections (pre-flight tests, noise, etc.), and trims the result to the
portion of the log where altitude is moving.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional, Tuple


def find_headers_and_peak(lines: List[str]) -> Tuple[int, Optional[int]]:
    """
    Returns (start_idx, end_idx) for the block that contains the max altitude.

    start_idx -> index of the header preceding the peak altitude
    end_idx   -> index of the header immediately after that block (or None)
    """
    if not lines:
        raise ValueError("Input file is empty.")

    canonical_header = lines[0].strip()
    header_indices: List[int] = []

    current_header_idx: Optional[int] = None
    peak_alt = float("-inf")
    peak_header_idx: Optional[int] = None
    peak_line_idx: Optional[int] = None

    for idx, raw in enumerate(lines):
        stripped = raw.strip()
        if not stripped:
            continue

        if stripped == canonical_header:
            header_indices.append(idx)
            current_header_idx = idx
            continue

        if current_header_idx is None:
            # Found data before any header; skip cautiously
            continue

        parts = stripped.split(",")
        if len(parts) <= 4:
            continue

        try:
            alt = float(parts[4])
        except ValueError:
            continue

        if alt > peak_alt:
            peak_alt = alt
            peak_header_idx = current_header_idx
            peak_line_idx = idx

    if peak_header_idx is None or peak_line_idx is None:
        raise ValueError("Could not determine peak altitude within the file.")

    next_header_idx = None
    for idx in header_indices:
        if idx > peak_line_idx:
            next_header_idx = idx
            break

    return peak_header_idx, next_header_idx


def extract_segment(lines: List[str], start_idx: int, end_idx: Optional[int]) -> List[str]:
    """
    Slice lines between start_idx and end_idx (or EOF if end_idx is None).
    """
    if end_idx is None:
        return lines[start_idx:]
    return lines[start_idx:end_idx]


def trim_to_motion(
    block_lines: List[str],
    diff_threshold: float,
    padding: int,
) -> List[str]:
    """
    Within a block (header + rows), keep only the portion where altitude changes.

    Args:
        block_lines: list containing exactly one header plus data rows
        diff_threshold: minimum absolute delta between consecutive altitude samples
                        that counts as "motion"
        padding: number of data rows to retain before/after detected motion
    """
    if len(block_lines) <= 1:
        return block_lines

    header = block_lines[0]
    data_lines = block_lines[1:]

    valid_records: List[Tuple[int, float]] = []
    for idx, raw in enumerate(data_lines):
        stripped = raw.strip()
        if not stripped:
            continue
        parts = stripped.split(",")
        if len(parts) <= 4:
            continue
        try:
            alt = float(parts[4])
        except ValueError:
            continue
        valid_records.append((idx, alt))

    if len(valid_records) < 2:
        return block_lines

    motion_start_idx: Optional[int] = None
    motion_end_idx: Optional[int] = None
    prev_alt = None
    prev_idx = None

    for idx, alt in valid_records:
        if prev_alt is not None and prev_idx is not None:
            if abs(alt - prev_alt) >= diff_threshold:
                if motion_start_idx is None:
                    motion_start_idx = prev_idx
                motion_end_idx = idx
        prev_alt = alt
        prev_idx = idx

    if motion_start_idx is None or motion_end_idx is None:
        # No meaningful motion detected; return block intact.
        return block_lines

    motion_start_idx = max(0, motion_start_idx - padding)
    motion_end_idx = min(len(data_lines) - 1, motion_end_idx + padding)

    trimmed_data = data_lines[motion_start_idx : motion_end_idx + 1]
    return [header] + trimmed_data


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Keep only the flight segment surrounding the highest altitude "
            "and trim to the timeframe where altitude is moving."
        )
    )
    parser.add_argument("input_csv", type=Path, help="Raw CSV containing multiple segments")
    parser.add_argument("output_csv", type=Path, help="Cleaned CSV output path")
    parser.add_argument(
        "--diff-threshold",
        type=float,
        default=0.2,
        help="Minimum |ΔAlt| between consecutive samples to count as motion (default: 0.2 m)",
    )
    parser.add_argument(
        "--padding",
        type=int,
        default=10,
        help="Number of samples to keep before/after detected motion (default: 10)",
    )
    args = parser.parse_args()

    if not args.input_csv.exists():
        raise FileNotFoundError(args.input_csv)

    with args.input_csv.open("r", encoding="utf-8") as fh:
        lines = fh.readlines()

    start_idx, end_idx = find_headers_and_peak(lines)
    block = extract_segment(lines, start_idx, end_idx)
    cleaned = trim_to_motion(block, args.diff_threshold, args.padding)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", encoding="utf-8") as fh:
        fh.writelines(cleaned)

    print(f"Identified peak segment: lines {start_idx} to {end_idx or len(lines)}")
    print(f"Saved cleaned data to {args.output_csv}")


if __name__ == "__main__":
    main()

