from constants import finalized_frame_list, not_main_finalized_frame_list, GOES_scale_value
from bisect import bisect_left
import numpy as np


def frame_sort_func(temp_frame_list, query_exposure_list, initial_value):
    total_index = 0
    print("temp frame list")
    print(temp_frame_list)
    for sublist in range(0, len(temp_frame_list)):
        sublist_length = len(temp_frame_list[sublist])
        for i in range(0, sublist_length):
            index_in_exposure_list = total_index + i + initial_value
            if i == sublist_length - 1:
                finalized_frame_list.append(index_in_exposure_list)
                total_index += sublist_length
                break
            if query_exposure_list[index_in_exposure_list] < 1.0:
                finalized_frame_list.append(index_in_exposure_list)
                total_index += sublist_length
                break



def frame_sort_func2(temp_frame_list, query_exposure_list, initial_value):
    total_index = 0
    print("temp frame list")
    print(temp_frame_list)
    for sublist in range(0, len(temp_frame_list)):
        sublist_length = len(temp_frame_list[sublist])
        for i in range(0, sublist_length):
            index_in_exposure_list = total_index + i + initial_value
            if i == sublist_length - 1:
                not_main_finalized_frame_list.append(index_in_exposure_list)
                total_index += sublist_length
                break
            if query_exposure_list[index_in_exposure_list] < 1.0:
                not_main_finalized_frame_list.append(index_in_exposure_list)
                total_index += sublist_length
                break

def split_frames(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


def take_closest(myList, myNumber):
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after
    else:
        return before


def findpeaks(series, DELTA):
    # Set inital values
    mn, mx = np.Inf, -np.Inf
    minpeaks = []
    maxpeaks = []
    lookformax = True
    start = True
    # Iterate over items in series
    for time_pos, value in series.iteritems():
        if value > mx:
            mx = value
            mxpos = time_pos
        if value < mn:
            mn = value
            mnpos = time_pos
        if lookformax:
            if value < mx - DELTA:
                # a local maxima
                maxpeaks.append((mxpos, mx * GOES_scale_value))
                mn = value
                mnpos = time_pos
                lookformax = False
            elif start:
                # a local minima at beginning
                minpeaks.append((mnpos, mn * GOES_scale_value))
                mx = value
                mxpos = time_pos
                start = False
        else:
            if value > mn + DELTA:
                # a local minima
                minpeaks.append((mnpos, mn * GOES_scale_value))
                mx = value
                mxpos = time_pos
                lookformax = True
    # check for extrema at end
    if value > mn + DELTA:
        maxpeaks.append((mxpos, mx * GOES_scale_value))
    elif value < mx - DELTA:
        minpeaks.append((mnpos, mn * GOES_scale_value))
    peak_time = max(maxpeaks)

    peak_time = str(max(maxpeaks))[12:-36].replace(" ", 'T')

    return minpeaks, maxpeaks