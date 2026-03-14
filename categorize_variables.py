"""
categorize_variables.py
-----------------------
Utilities for categorizing variables by data type and statistical role.
"""

from __future__ import annotations

import math
from typing import Any


# ---------------------------------------------------------------------------
# Type-based categorization
# ---------------------------------------------------------------------------

def categorize_by_type(variable: Any) -> str:
    """Return the category name that best describes *variable*'s Python type.

    Categories
    ----------
    - "integer"   : int (but not bool)
    - "float"     : float
    - "complex"   : complex
    - "boolean"   : bool
    - "string"    : str
    - "bytes"     : bytes / bytearray
    - "sequence"  : list / tuple
    - "mapping"   : dict
    - "set"       : set / frozenset
    - "none"      : None
    - "unknown"   : anything else

    Examples
    --------
    >>> categorize_by_type(42)
    'integer'
    >>> categorize_by_type(3.14)
    'float'
    >>> categorize_by_type("hello")
    'string'
    >>> categorize_by_type(True)
    'boolean'
    >>> categorize_by_type([1, 2, 3])
    'sequence'
    >>> categorize_by_type(None)
    'none'
    """
    # bool is a subclass of int – check it first
    if isinstance(variable, bool):
        return "boolean"
    if isinstance(variable, int):
        return "integer"
    if isinstance(variable, float):
        return "float"
    if isinstance(variable, complex):
        return "complex"
    if isinstance(variable, str):
        return "string"
    if isinstance(variable, (bytes, bytearray)):
        return "bytes"
    if isinstance(variable, (list, tuple)):
        return "sequence"
    if isinstance(variable, dict):
        return "mapping"
    if isinstance(variable, (set, frozenset)):
        return "set"
    if variable is None:
        return "none"
    return "unknown"


# ---------------------------------------------------------------------------
# Statistical / measurement-scale categorization
# ---------------------------------------------------------------------------

def categorize_statistical(values: list[Any]) -> str:
    """Infer the statistical measurement scale of a list of values.

    Scales (in order of richness)
    ------------------------------
    - "nominal"  : discrete labels with no natural order (strings, mixed types)
    - "ordinal"  : discrete labels that can be sorted (integers/floats with
                   fewer than ORDINAL_THRESHOLD distinct values)
    - "interval" : numeric, many distinct values, no meaningful zero
                   (not detected automatically – returned as "ratio")
    - "ratio"    : numeric, many distinct values (default for numeric data)
    - "binary"   : exactly two distinct values
    - "empty"    : the list is empty

    Parameters
    ----------
    values:
        A non-empty list of raw values.

    Examples
    --------
    >>> categorize_statistical([1, 0, 1, 0, 1])
    'binary'
    >>> categorize_statistical(["cat", "dog", "bird"])
    'nominal'
    >>> categorize_statistical([1, 2, 3, 1, 2])
    'ordinal'
    >>> categorize_statistical([1.5, 2.7, 3.14, 100.0, 0.001, 42.0])
    'ratio'
    """
    ORDINAL_THRESHOLD = 20  # max distinct values to still call it ordinal

    if not values:
        return "empty"

    distinct = set(values)

    if len(distinct) == 2:
        return "binary"

    # Determine whether values are all numeric
    all_numeric = all(isinstance(v, (int, float)) and not isinstance(v, bool)
                      for v in values)

    if not all_numeric:
        return "nominal"

    if len(distinct) <= ORDINAL_THRESHOLD:
        return "ordinal"

    return "ratio"


# ---------------------------------------------------------------------------
# Numeric range categorization
# ---------------------------------------------------------------------------

def categorize_numeric_range(value: int | float) -> str:
    """Classify a numeric value into a named range bucket.

    Buckets
    -------
    - "nan"       : not a number (float NaN)
    - "negative"  : value < 0
    - "zero"      : value == 0
    - "small"     : 0 < value <= 1
    - "medium"    : 1 < value <= 1_000
    - "large"     : 1_000 < value <= 1_000_000
    - "very_large": value > 1_000_000
    - "infinite"  : positive or negative infinity

    Examples
    --------
    >>> categorize_numeric_range(0)
    'zero'
    >>> categorize_numeric_range(-5)
    'negative'
    >>> categorize_numeric_range(0.5)
    'small'
    >>> categorize_numeric_range(500)
    'medium'
    >>> categorize_numeric_range(1_500_000)
    'very_large'
    """
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        if math.isinf(value):
            return "infinite"
    if value < 0:
        return "negative"
    if value == 0:
        return "zero"
    if value <= 1:
        return "small"
    if value <= 1_000:
        return "medium"
    if value <= 1_000_000:
        return "large"
    return "very_large"


# ---------------------------------------------------------------------------
# Convenience: categorize a whole dictionary of named variables
# ---------------------------------------------------------------------------

def categorize_dict(variables: dict[str, Any]) -> dict[str, dict[str, str]]:
    """Categorize every value in a name→value dictionary.

    Returns a dict mapping each name to a sub-dict with keys:
    - ``"type_category"``  from :func:`categorize_by_type`
    - ``"range_category"`` from :func:`categorize_numeric_range`
                           (only for numeric scalars, else ``"n/a"``)

    Examples
    --------
    >>> result = categorize_dict({"age": 25, "name": "Alice", "score": 98.6})
    >>> result["age"]
    {'type_category': 'integer', 'range_category': 'medium'}
    >>> result["name"]
    {'type_category': 'string', 'range_category': 'n/a'}
    """
    output: dict[str, dict[str, str]] = {}
    for name, value in variables.items():
        type_cat = categorize_by_type(value)
        if type_cat in ("integer", "float") and not isinstance(value, bool):
            range_cat = categorize_numeric_range(value)
        else:
            range_cat = "n/a"
        output[name] = {
            "type_category": type_cat,
            "range_category": range_cat,
        }
    return output


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    sample_variables = {
        "age":        25,
        "height":     1.75,
        "name":       "Alice",
        "active":     True,
        "score":      float("nan"),
        "distance":   1_500_000.0,
        "tags":       ["python", "data"],
        "config":     {"debug": False},
        "nothing":    None,
    }

    print("=== Type & Range Categorization ===")
    results = categorize_dict(sample_variables)
    for var_name, cats in results.items():
        print(f"  {var_name:<12} -> {cats}")

    print()
    print("=== Statistical Scale (example lists) ===")
    examples = {
        "binary flags":   [0, 1, 1, 0, 1],
        "grade levels":   [1, 2, 3, 2, 1, 3],
        "temperatures":   [36.6, 37.1, 36.9, 38.0, 37.5, 36.8, 39.2],
        "city names":     ["Rome", "Milan", "Naples", "Turin"],
        "empty list":     [],
    }
    for label, vals in examples.items():
        print(f"  {label:<20} -> {categorize_statistical(vals)}")
