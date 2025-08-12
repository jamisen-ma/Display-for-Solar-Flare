from __future__ import annotations

import bisect
from typing import List, Sequence

from PyQt5 import QtCore

from logging_setup import get_logger


logger = get_logger(__name__)


def normalized_plot(plot_list: Sequence[float]) -> List[float]:
    max_value = max(plot_list) if plot_list else 1.0
    return [float(intensity) / max_value for intensity in plot_list]


def selected_RHESSI_keys(flare_list: Sequence[str], selected_start_date: str, selected_end_date: str) -> List[str]:
    left_index = bisect.bisect(flare_list, selected_start_date) - 1
    right_index = bisect.bisect(flare_list, selected_end_date)
    return list(flare_list[left_index:right_index])


def clickBox(self, state):
    logger.debug('Checkbox state: %s', 'Checked' if state == QtCore.Qt.Checked else 'Unchecked')