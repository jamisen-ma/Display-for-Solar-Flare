import bisect
from PyQt5 import QtCore

def normalized_plot(plot_list):
    new_list = []
    max_value = max(plot_list)
    for intensity in plot_list:
        new_val = intensity / max_value
        new_list.append(new_val)
    return new_list


def selected_RHESSI_keys(flare_list, selected_start_date, selected_end_date):
    left_index = bisect.bisect(flare_list, selected_start_date) - 1

    right_index = bisect.bisect(flare_list, selected_end_date)

    selected_solar_flares_keys_list = flare_list[left_index:right_index]

    return selected_solar_flares_keys_list


def clickBox(self, state):
    if state == QtCore.Qt.Checked:
        print('Checked')
    else:
        print('Unchecked')