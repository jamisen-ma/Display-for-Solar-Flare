import pandas as pd
from datetime import timedelta
import numpy as np
import bisect
from units_conversion_functions import *
from math_time_functions import *
from manipulate_flare_lists_functions import frame_sort_func, split_frames, take_closest
from HMI_data_analysis import get_HMI_data, get_pixel_coords_from_AR_center
from PyQt_shapes import RectItem, EclipseItem
from image_processing_functions import strided_rescale
from creating_image_lists_functions import creating_1600_image_list, create_difference_and_ratio_lists, create_date_list
from get_online_files import get_image_urls

smaller_solar_dict = {}


date_format2 = "%Y-%m-%d %H:%M:%S"

dates_for_searching = create_date_list(start_date, end_date)
url_list = get_image_urls(dates_for_searching)

solar_flares_panda_frame = pd.read_csv('hessi.solar.flare.UP_To_2018.csv')
# solar_flares_panda_frame['image_data'] = np.nan

solar_flares_panda_frame['flag.2'] = solar_flares_panda_frame['flag.2'].fillna("")
solar_flares_panda_frame['flag.3'] = solar_flares_panda_frame['flag.3'].fillna("")
solar_flares_panda_frame['flag.4'] = solar_flares_panda_frame['flag.4'].fillna("")
solar_flares_panda_frame['flag.5'] = solar_flares_panda_frame['flag.5'].fillna("")

solar_flares_panda_frame['start_time_placeholder'] = solar_flares_panda_frame['start.time']
solar_flares_panda_frame['start.time'] = pd.to_datetime(solar_flares_panda_frame['start.time'], format='%H:%M:%S')
solar_flares_panda_frame['start.date'] = pd.to_datetime(solar_flares_panda_frame['start.date'], format="%Y-%m-%d")

solar_flares_panda_frame['peak_placeholder'] = solar_flares_panda_frame['peak']
solar_flares_panda_frame['peak'] = pd.to_datetime(solar_flares_panda_frame['peak'], format='%H:%M:%S')
solar_flares_panda_frame['real.start.date'] = np.where(
    (solar_flares_panda_frame['start.time'] > solar_flares_panda_frame['peak']),
    solar_flares_panda_frame['start.date'] + timedelta(days=1), solar_flares_panda_frame['start.date'])

solar_flares_panda_frame['real.end.date'] = np.where(
    (solar_flares_panda_frame['end'] < solar_flares_panda_frame['start_time_placeholder']),
    solar_flares_panda_frame['start.date'] + timedelta(days=1), solar_flares_panda_frame['start.date'])

solar_flares_panda_frame['real.peak.date'] = np.where(
    (solar_flares_panda_frame['peak_placeholder'] < solar_flares_panda_frame['start_time_placeholder']),
    solar_flares_panda_frame['start.date'] + timedelta(days=1), solar_flares_panda_frame['start.date'])

solar_flares_panda_frame['real_start_time'] = solar_flares_panda_frame['real.start.date'].dt.strftime(
    '%Y-%m-%d') + " " + solar_flares_panda_frame['start_time_placeholder']
print("finsiehedddd")
solar_flares_panda_frame['real_peak_time'] = solar_flares_panda_frame['real.peak.date'].dt.strftime(
    '%Y-%m-%d') + " " + solar_flares_panda_frame['peak_placeholder']

solar_flares_panda_frame['real_end_time'] = solar_flares_panda_frame['real.end.date'].dt.strftime(
    '%Y-%m-%d') + " " + \
                                            solar_flares_panda_frame['end']
for x in solar_flares_panda_frame.select_dtypes(include=['datetime64']).columns.tolist():
    solar_flares_panda_frame[x] = solar_flares_panda_frame[x].astype(str)

solar_flares_panda_frame['real.start.date'] = solar_flares_panda_frame['real.start.date']
solar_flares_panda_frame['peak'] = solar_flares_panda_frame['peak'].str[11:]
solar_flares_panda_frame["peak_time"] = solar_flares_panda_frame['real.start.date'] + " " + \
                                        solar_flares_panda_frame['peak']
solar_flares_panda_frame["flags"] = solar_flares_panda_frame['flag.1'].astype(str) + solar_flares_panda_frame[
    'flag.2'] + solar_flares_panda_frame['flag.3'] + solar_flares_panda_frame['flag.4'] + solar_flares_panda_frame[
                                        'flag.5']

solar_flares_panda_frame = solar_flares_panda_frame.rename(
    columns={'x.pos.asec': 'RHESSI List X center', 'y.pos.asec': 'RHESSI List Y center',
             'real_start_time': 'RHESSI Start Time', 'real_peak_time': 'RHESSI Peak Time',
             'real_end_time': 'RHESSI End Time'})

start_time = solar_flares_panda_frame['RHESSI Start Time']
solar_flares_panda_frame.drop(labels=['RHESSI Start Time'], axis=1, inplace=True)
solar_flares_panda_frame.insert(0, 'RHESSI Start Time', start_time)
peak_time = solar_flares_panda_frame['RHESSI Peak Time']
solar_flares_panda_frame.drop(labels=['RHESSI Peak Time'], axis=1, inplace=True)
solar_flares_panda_frame.insert(1, 'RHESSI Peak Time', peak_time)
end_time = solar_flares_panda_frame['RHESSI End Time']
solar_flares_panda_frame.drop(labels=['RHESSI End Time'], axis=1, inplace=True)
solar_flares_panda_frame.insert(2, 'RHESSI End Time', end_time)


def make_new_column_func(key, big_list, data_frame):
    data_frame[key] = data_frame[big_list].str.contains(key, na=True)


flag_list = ['a0', 'a1', 'a2', 'a3', 'A0', 'A1', 'A2', 'A3', 'DF', 'DR', 'ED', 'EE', 'ES', 'FE', 'FR', 'FS', 'GD',
             'GE', 'MR', 'NS', 'PE', 'PS', 'SD', 'SE', 'SS']

solar_flares_panda_frame = solar_flares_panda_frame.drop(
    ['real.start.date', 'real.end.date', 'real.peak.date', 'radial', 'start.date', 'peak', 'end', 'start.time'],
    axis=1)
print("Start time list")
# solar_flares_panda_frame = solar_flares_panda_frame.drop(solar_flares_panda_frame.index[range(32,116143)])
print(solar_flares_panda_frame.index)
print(start_time_list)
converted_to_string_start_time_list = []

for date2 in start_time_list:
    # date2 = datetime.strptime(str(date2), '%Y-%m-%d %H:%M:%S')
    date2 = str(date2)
    converted_to_string_start_time_list.append(date2)
print(converted_to_string_start_time_list)
# solar_flares_panda_frame['GOES Start Time'] = [1,2,3,4,5,6,7]
for flag in flag_list:
    make_new_column_func(flag, 'flags', solar_flares_panda_frame)

solar_flares_dict = solar_flares_panda_frame.set_index("peak_time").T.to_dict('dict')


def selected_RHESSI_keys(flare_list, selected_start_date, selected_end_date):
    left_index = bisect.bisect(flare_list, selected_start_date)
    right_index = bisect.bisect(flare_list, selected_end_date)
    print("numbers")
    print(left_index, right_index)
    selected_solar_flares_keys_list = flare_list[left_index:right_index]

    return selected_solar_flares_keys_list


solar_flares_keys_list = list(solar_flares_dict.keys())

solar_flares_time_list = selected_RHESSI_keys(solar_flares_keys_list, start_date, end_date)

for i in range(0, len(solar_flares_time_list)):
    flare_info = solar_flares_dict[solar_flares_time_list[i]]
    smaller_solar_dict[solar_flares_time_list[i]] = flare_info

print("smaller_solar_dict")
print(solar_flares_time_list)


# ________

def match_RHESSI_image_with_RHESSI_list(flare_list, urls):
    print("urls")
    print(urls)
    for y in range(0, len(urls)):

        print(y)

        headers = fits.open(urls[y])[0].header

        peak_time = headers['DATE_OBS'].replace("T", " ")

        flare_index = bisect.bisect_left(flare_list, peak_time)

        peak_time = to_datetime_object(str(peak_time)[:-4], date_format2)
        print(peak_time)
        print(flare_index)

        try:
            if abs(to_datetime_object(flare_list[flare_index], date_format2) - peak_time) > abs(
                    to_datetime_object(flare_list[flare_index + 1], date_format2) - peak_time):
                flare_dict_key = flare_list[flare_index + 1]

            else:
                flare_dict_key = flare_list[flare_index]
        except IndexError:
            continue

        # image_data = fits.getdata(urls[y], ext=0)
        print("image received")
        rhessi_header_x = headers['CRPIX1']
        rhessi_header_y = headers['CRPIX1']

        smaller_solar_dict[flare_dict_key]['rhessi_header_x'] = rhessi_header_x
        smaller_solar_dict[flare_dict_key]['rhessi_header_y'] = rhessi_header_y
        smaller_solar_dict[flare_dict_key]['rhessi_conversion'] = headers['CDELT1']
        # smaller_solar_dict[flare_dict_key]['image_data'] = image_data
        smaller_solar_dict[flare_dict_key]['rhessi_image_link'] = urls[y]
        print("flare dict key")
        print(flare_dict_key)


match_RHESSI_image_with_RHESSI_list(solar_flares_time_list, url_list)
solar_flares_datetime_object_time_list = [to_datetime_object(date_val, '%Y-%m-%d %H:%M:%S') for date_val in
                                          solar_flares_time_list]

print("dictssss")
print(solar_flares_datetime_object_time_list)
print(real_peak_time_list)
print("real time list")
print(real_start_time_list)
for i in range(0, len(real_peak_time_list)):
    flare_index = bisect.bisect(solar_flares_datetime_object_time_list, real_peak_time_list[i])

    try:
        if abs(solar_flares_datetime_object_time_list[flare_index - 1] - real_peak_time_list[i]) < timedelta(
                minutes=20):
            data = smaller_solar_dict[solar_flares_time_list[flare_index - 1]]
            master_dict[str(real_peak_time_list[i])] = data

            master_dict[str(real_peak_time_list[i])]['flare.class'] = real_flare_class_list[i]
            master_dict[str(real_peak_time_list[i])]['GOES Start Time'] = real_start_time_list[i]
            master_dict[str(real_peak_time_list[i])]['GOES Peak Time'] = real_peak_time_list[i]
            master_dict[str(real_peak_time_list[i])]['GOES End Time'] = real_end_time_list[i]
            print("success")

        elif abs(solar_flares_datetime_object_time_list[flare_index] - real_peak_time_list[i]) < timedelta(
                minutes=20):
            print("success1")
            data = smaller_solar_dict[solar_flares_time_list[flare_index]]
            master_dict[str(real_peak_time_list[i])] = data
            master_dict[str(real_peak_time_list[i])]['flare.class'] = real_flare_class_list[i]
            master_dict[str(real_peak_time_list[i])]['GOES Start Time'] = real_start_time_list[i]
            master_dict[str(real_peak_time_list[i])]['GOES Peak Time'] = real_peak_time_list[i]
            master_dict[str(real_peak_time_list[i])]['GOES End Time'] = real_end_time_list[i]

    except IndexError:
        pass
print("master dict 2")

new_dataframe = pd.DataFrame(master_dict).transpose()
drop_list = ['flare', 'flags', 'a0', 'a1', 'a2', 'a3',
             'A0', 'A1', 'A2', 'A3', 'DF', 'DR', 'ED', 'EE', 'ES', 'FE', 'FR', 'FS',
             'GD', 'GE', 'MR', 'NS', 'PE', 'PS', 'SD', 'SE', 'SS', 'flag.1', 'flag.2', 'flag.3', 'flag.4',
             'flag.5', 'energy.kev', 'start_time_placeholder', 'peak_placeholder']

for element in drop_list:
    new_dataframe.drop(element, axis=1, inplace=True)
print(new_dataframe)
print(new_dataframe.index)
print(new_dataframe.columns)
flare_class1 = new_dataframe['flare.class']
new_dataframe.drop(labels=['flare.class'], axis=1, inplace=True)
new_dataframe.insert(3, 'Flare Class', flare_class1)
pd.to_pickle(new_dataframe, 'dataframetest.pickle')

print('real peak time list')
print(real_peak_time_list)