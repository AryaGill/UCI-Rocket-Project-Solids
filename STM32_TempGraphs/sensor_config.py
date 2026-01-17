"""
Sensor configuration file - add new sensors here
"""

SENSORS = {
    'LPS22HHTR': {
        'columns': ['pressure', 'temperature', 'altitude'],
        'graphs': [
            {
                'title': 'Pressure',
                'columns': ['pressure'],
                'y_label': 'Pressure (hPa)',
                'colors': ['#FF6B6B']  # Red
            },
            {
                'title': 'Temperature',
                'columns': ['temperature'],
                'y_label': 'Temperature (°C)',
                'colors': ['#4ECDC4']  # Cyan
            },
            {
                'title': 'Altitude',
                'columns': ['altitude'],
                'y_label': 'Altitude (m)',
                'colors': ['#45B7D1']  # Blue
            }
        ],
        'window_size': 100  # Number of data points to display
    }
}

# Color mapping for consistent colors across all graphs
AXIS_COLORS = {
    'pressure': '#FF6B6B',     # Red
    'temperature': '#4ECDC4',  # Cyan
    'altitude': '#45B7D1'      # Blue
}