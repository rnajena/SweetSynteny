# (C) 2024, Maria Schreiber, MIT license
"""Reading color.json"""
import json
import os

DEFAULT_COLORS = {
    "goi": "#E63946",
    "bars": "#457B9D",
    "lines": "#457B9D",
    "histogram": "#457B9D",
    "scatter": "#457B9D",
    "heatmap_absent": "#FFFFFF",
    "heatmap_present": "#333333",
    "background": "#FFFFFF",
    "text": "#222222",
    "grid": "#D9D9D9",
}


def load_colors(filename=None):
    """Load a user palette, filling unspecified keys with defaults."""
    colors = DEFAULT_COLORS.copy()

    if filename:
        if not os.path.isfile(filename):
            raise FileNotFoundError(
                f"Color scheme file not found: {filename}"
            )
        with open(filename, encoding="utf-8") as handle:
            supplied = json.load(handle)

        if not isinstance(supplied, dict):
            raise ValueError("Color scheme must be a JSON object.")

        unknown = set(supplied) - set(colors)
        if unknown:
            raise ValueError(
                f"Unknown color keys: {sorted(unknown)}"
            )

        colors.update(supplied)

    return colors
