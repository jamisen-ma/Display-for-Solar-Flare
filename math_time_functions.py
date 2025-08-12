"""Math and time utility helpers."""

from __future__ import annotations

import math
from datetime import datetime, timedelta
from typing import Any


def to_datetime_object(date_string: Any, date_format: str) -> datetime:
    """Parse a date string to a ``datetime`` using the provided format.

    Accepts any value that can be stringified to the expected format.
    """
    return datetime.strptime(str(date_string), date_format)


def convert_timedelta(duration: timedelta) -> int:
    """Convert ``timedelta`` to rounded-up minutes (minimum step of 1 minute)."""
    seconds = int(duration.total_seconds())
    return max(1, math.ceil(seconds / 60))


def roundup(x: float) -> int:
    """Round x up to the nearest 100."""
    return int(math.ceil(x / 100.0)) * 100


def rounddown(x: float) -> int:
    """Round x down to the nearest 100."""
    return (int(math.ceil(x / 100.0)) * 100) - 100