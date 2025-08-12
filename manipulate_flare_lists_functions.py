from __future__ import annotations

from typing import Generator, Iterable, List, Sequence, Tuple

from bisect import bisect_left
import numpy as np

from constants import finalized_frame_list, not_main_finalized_frame_list, GOES_scale_value


def frame_sort_func(temp_frame_list: Sequence[Sequence[int]], query_exposure_list: Sequence[float], initial_value: int) -> None:
    total_index = 0
    for sublist in range(0, len(temp_frame_list)):
        sublist_length = len(temp_frame_list[sublist])
        for i in range(0, sublist_length):
            index_in_exposure_list = total_index + i + initial_value
            if i == sublist_length - 1 or query_exposure_list[index_in_exposure_list] < 1.0:
                finalized_frame_list.append(index_in_exposure_list)
                total_index += sublist_length
                break


def frame_sort_func2(temp_frame_list: Sequence[Sequence[int]], query_exposure_list: Sequence[float], initial_value: int) -> None:
    total_index = 0
    for sublist in range(0, len(temp_frame_list)):
        sublist_length = len(temp_frame_list[sublist])
        for i in range(0, sublist_length):
            index_in_exposure_list = total_index + i + initial_value
            if i == sublist_length - 1 or query_exposure_list[index_in_exposure_list] < 1.0:
                not_main_finalized_frame_list.append(index_in_exposure_list)
                total_index += sublist_length
                break


def split_frames(a: Sequence[int], n: int) -> Generator[Sequence[int], None, None]:
    k, m = divmod(len(a), n)
    for i in range(n):
        yield a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)]


def take_closest(myList: Sequence[float], myNumber: float) -> float:
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    return after if after - myNumber < myNumber - before else before


def findpeaks(series, DELTA: float):
    mn, mx = np.Inf, -np.Inf
    minpeaks = []
    maxpeaks = []
    lookformax = True
    start = True
    for time_pos, value in series.iteritems():
        if value > mx:
            mx = value
            mxpos = time_pos
        if value < mn:
            mn = value
            mnpos = time_pos
        if lookformax:
            if value < mx - DELTA:
                maxpeaks.append((mxpos, mx * GOES_scale_value))
                mn = value
                mnpos = time_pos
                lookformax = False
            elif start:
                minpeaks.append((mnpos, mn * GOES_scale_value))
                mx = value
                mxpos = time_pos
                start = False
        else:
            if value > mn + DELTA:
                minpeaks.append((mnpos, mn * GOES_scale_value))
                mx = value
                mxpos = time_pos
                lookformax = True
    if value > mn + DELTA:
        maxpeaks.append((mxpos, mx * GOES_scale_value))
    elif value < mx - DELTA:
        minpeaks.append((mnpos, mn * GOES_scale_value))
    return minpeaks, maxpeaks