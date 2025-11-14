#!/usr/bin/env python3
"""
Simple CSV plotting utility.

Reads a CSV file, handles duplicate headers, and plots specified columns.

Usage:
    python basic_plot.py path/to/file.csv Alt Time
    python basic_plot.py path/to/file.csv Alt --x-axis Time
    python basic_plot.py path/to/file.csv Alt Temp Press --x-axis Time
"""
import argparse
from pathlib import Path
from typing import List, Optional

import matplotlib.pyplot as plt
import pandas as pd


def read_csv_with_duplicate_headers(csv_path: Path) -> pd.DataFrame:
    """
    Read CSV file, skipping duplicate header rows.
    
    Detects rows that match the first header and skips them.
    """
    # First, read just to get the header
    with open(csv_path, 'r') as f:
        first_line = f.readline().strip()
        header = first_line.split(',')
    
    # Read the file and filter out duplicate headers
    rows = []
    with open(csv_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            stripped = line.strip()
            # Skip empty lines
            if not stripped:
                continue
            # Skip lines that match the header exactly
            if stripped == first_line:
                if line_num == 1:
                    # Keep the first header
                    rows.append(line)
                # Skip subsequent duplicate headers
                continue
            rows.append(line)
    
    # Parse the filtered rows
    from io import StringIO
    csv_content = ''.join(rows)
    df = pd.read_csv(StringIO(csv_content))
    
    return df


def plot_columns(
    df: pd.DataFrame,
    columns: List[str],
    x_axis: Optional[str] = None,
    output_path: Optional[Path] = None,
) -> None:
    """
    Plot specified columns from the dataframe.
    
    Args:
        df: DataFrame containing the data
        columns: List of column names to plot
        x_axis: Column name to use for x-axis (default: index or 'Time' if available)
        output_path: Optional path to save the plot
    """
    # Determine x-axis
    if x_axis:
        if x_axis not in df.columns:
            raise ValueError(f"X-axis column '{x_axis}' not found in CSV. Available columns: {list(df.columns)}")
        x_data = df[x_axis]
        x_label = x_axis
        # Convert time from ms to seconds if it looks like milliseconds
        if x_axis.lower() == 'time' and x_data.max() > 10000:
            x_data = x_data / 1000.0
            x_label = "Time (s)"
    else:
        # Try to use Time column if available
        if 'Time' in df.columns:
            x_data = df['Time']
            if x_data.max() > 10000:
                x_data = x_data / 1000.0
            x_label = "Time (s)"
        else:
            x_data = df.index
            x_label = "Sample Index"
    
    # Check all requested columns exist
    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise ValueError(
            f"Columns not found: {missing}. Available columns: {list(df.columns)}"
        )
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    for col in columns:
        ax.plot(x_data, df[col], label=col, marker='.', markersize=1, alpha=0.7)
    
    ax.set_xlabel(x_label)
    ax.set_ylabel("Value")
    ax.set_title(f"CSV Data: {', '.join(columns)}")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300)
        print(f"Plot saved to: {output_path}")
    else:
        plt.show()
    
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot columns from a CSV file, handling duplicate headers."
    )
    parser.add_argument(
        "csv_path",
        type=Path,
        help="Path to the CSV file"
    )
    parser.add_argument(
        "columns",
        nargs="+",
        help="Column names to plot (one or more)"
    )
    parser.add_argument(
        "--x-axis",
        type=str,
        default=None,
        help="Column to use for x-axis (default: 'Time' if available, otherwise index)"
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=None,
        help="Output path for the plot (PNG). If not specified, displays interactively."
    )
    args = parser.parse_args()
    
    if not args.csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {args.csv_path}")
    
    print(f"Reading CSV: {args.csv_path}")
    df = read_csv_with_duplicate_headers(args.csv_path)
    print(f"Loaded {len(df)} rows, {len(df.columns)} columns")
    print(f"Columns: {', '.join(df.columns)}")
    
    plot_columns(df, args.columns, x_axis=args.x_axis, output_path=args.output)


if __name__ == "__main__":
    main()

