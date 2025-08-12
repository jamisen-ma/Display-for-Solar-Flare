from __future__ import absolute_import, division, print_function
import matplotlib
import scipy.ndimage

matplotlib.use("TkAgg")
from matplotlib import style
import requests
from datetime import date, timedelta
import pandas as pd
import json
import codecs
import bisect
from tkinter import ttk
import datetime
from numpy.lib.stride_tricks import as_strided
import re
import matplotlib.backends.backend_tkagg as tkagg
import matplotlib.ticker as mticker
from astropy.time import Time

from matplotlib.backend_bases import Event
import matplotlib.ticker as mtick

import gc
import numpy as np
import sys
from sunpy.map import Map
import drms
from matplotlib.dates import YearLocator
import astropy.units as u

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.patches as patches
import math
from astropy.visualization import AsinhStretch
from matplotlib.animation import FuncAnimation
from astropy.visualization.mpl_normalize import ImageNormalize
import matplotlib.pyplot as plt
import sunpy.timeseries as ts
from astropy.io import fits
from astropy.time import TimeDelta
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net import hek
from sunpy.time import TimeRange, parse_time
from sunpy.timeseries import TimeSeries
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from tkinter import *
from PIL import ImageTk, Image
from PIL import Image
from astropy.table import Table
import sunpy.map
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
import matplotlib.lines as lines
import timeit
from bs4 import BeautifulSoup
import requests
import pandas as pd
import time as lmao

# Centralized logging and shared utilities
from logging_setup import get_logger
from image_processing_functions import (
    strided_rescale as _strided_rescale,
    mask_saturation as _mask_saturation,
)
from math_time_functions import (
    to_datetime_object as _to_datetime_object,
    convert_timedelta as _convert_timedelta,
)
from creating_image_lists_functions import (
    create_date_list as _create_date_list,
)
from get_online_files import (
    get_image_urls as _get_image_urls,
)
from units_conversion_functions import (
    x_helioprojective_to_pixel as _x_hp_to_px,
    y_helioprojective_to_pixel as _y_hp_to_px,
    x_pixel_to_helioprojective as _x_px_to_hp,
    y_pixel_to_helioprojective as _y_px_to_hp,
    rhessi_pixel_to_helioprojective_x as _r_px_to_hp_x,
    rhessi_pixel_to_helioprojective_y as _r_px_to_hp_y,
    rhessi_x_helioprojective_to_pixel as _r_hp_to_px_x,
    rhessi_y_helioprojective_to_pixel as _r_hp_to_px_y,
    helioprojective_to_heliographic as _hp_to_hg,
)

logger = get_logger(__name__)

space = 155

style.use("seaborn-bright")

root = Tk()
root.title("Finding Solar Flares")
root.geometry("400x350")
root.configure(background='white')

odd = True
finalized_frame_list = []
full_image_list_dict = []
sun_center_x_list_placeholder = []
sun_center_y_list_placeholder = []
pixel_arcsec_x_list_placeholder = []
pixel_arcsec_y_list_placeholder = []
auto_verification = []
AEC_list = []
key = []
manual_selected_time_range = []
GOES_graph_placeholder_list = ["logarithmic"]
wavelength_list = [94, 131, 171, 193, 211, 304, 335]
current_frame_number = 0
flare_selected = False
one_frame_state = False
graph_type_regular_bool = True
graph_type_changed_bool = False
GOES_frame = True
slider_changed = False
animation_interval_changed = False
user_clicked_map_bool = False
slider_used_func = False
coord_type_changed_bool = True
previous_flare_dict = []
real_peak_time_list = []
real_flare_class_list = []
flares_for_rhessi_comparison = []
flares_dropdown_menu = []
master_dict = {}
smaller_solar_dict = {}
rhessi_axis_dict = {}

animation_interval = 25
min_value = 0
max_value_vmax = 10000
max_threshold = 4000

label_distance = 20
initial_index = 0
GOES_scale_value = 0.99

cropped_dimension = 15

jsoc = drms.Client()
ds = jsoc.series(r"aia.lev1_euv_12s")  # Level 1.0
si = jsoc.info(ds[0])


def interactive_legend(ax=None):
    if ax is None:
        ax = plt.gca()
    if ax.legend_ is None:
        ax.legend()

    return InteractiveLegend(ax.get_legend())


class InteractiveLegend(object):
    def __init__(self, legend):
        self.legend = legend
        self.fig = legend.axes.figure

        self.lookup_artist, self.lookup_handle = self._build_lookups(legend)
        self._setup_connections()

        self.update()

    def _setup_connections(self):
        for artist in self.legend.texts + self.legend.legendHandles:
            artist.set_picker(10)  # 10 points tolerance

        self.fig.canvas.mpl_connect('pick_event', self.on_pick)
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)

    def _build_lookups(self, legend):
        labels = [t.get_text() for t in legend.texts]
        handles = legend.legendHandles
        label2handle = dict(zip(labels, handles))
        handle2text = dict(zip(handles, legend.texts))

        lookup_artist = {}
        lookup_handle = {}
        for artist in legend.axes.get_children():
            if artist.get_label() in labels:
                handle = label2handle[artist.get_label()]
                lookup_handle[artist] = handle
                lookup_artist[handle] = artist
                lookup_artist[handle2text[handle]] = artist

        lookup_handle.update(zip(handles, handles))
        lookup_handle.update(zip(legend.texts, handles))

        return lookup_artist, lookup_handle

    def on_pick(self, event):
        handle = event.artist
        if handle in self.lookup_artist:
            artist = self.lookup_artist[handle]
            artist.set_visible(not artist.get_visible())
            self.update()

    def on_click(self, event):
        if event.button == 3:
            visible = False
        elif event.button == 2:
            visible = True
        else:
            return

        for artist in self.lookup_artist.values():
            artist.set_visible(visible)
        self.update()

    def update(self):
        for artist in self.lookup_artist.values():
            handle = self.lookup_handle[artist]
            if artist.get_visible():
                handle.set_visible(True)
            else:
                handle.set_visible(False)
        self.fig.canvas.draw()

    def show(self):
        plt.show()


def do_nothing():
    pass


def create_frame():
    global frame1
    frame1 = Frame(root, width=670, height=900, bg="white", padx=0, pady=0, highlightthickness=0, borderwidth=0)
    frame1.grid(row=0, column=0)
    root.grid_propagate(False)


def strided_rescale(g, bin_fac):
    return _strided_rescale(g, bin_fac)


create_frame()


def create_frame2():
    global frame2
    frame2 = Frame(root, width=670, height=900, bg="white", padx=0, pady=0, highlightthickness=0, borderwidth=0)
    frame2.pack(side="top", anchor=NW)

    root.grid_propagate(False)


def roundup(x):
    return int(math.ceil(x / 100.0)) * 100


def rounddown(x):
    return (int(math.ceil(x / 100.0)) * 100) - 100


def x_helioprojective_to_pixel(value, additional_value):
    return _x_hp_to_px(value, additional_value, pixel_arcsec_x_list_placeholder, sun_center_x_list_placeholder)


def y_helioprojective_to_pixel(value, additional_value):
    return _y_hp_to_px(value, additional_value, pixel_arcsec_y_list_placeholder, sun_center_y_list_placeholder)


def x_pixel_to_helioprojective(value):
    return _x_px_to_hp(value, sun_center_x_list_placeholder, pixel_arcsec_x_list_placeholder)


def y_pixel_to_helioprojective(value):
    return _y_px_to_hp(value, sun_center_y_list_placeholder, pixel_arcsec_y_list_placeholder)


def rhessi_pixel_to_helioprojective_x(header_center, value, multiplier):
    return _r_px_to_hp_x(header_center, value, multiplier)


def rhessi_pixel_to_helioprojective_y(header_center, value, multiplier):
    return _r_px_to_hp_y(header_center, value, multiplier)


def rhessi_x_helioprojective_to_pixel(header_center, value, multiplier, additional_value):
    return _r_hp_to_px_x(header_center, value, multiplier, additional_value)


def rhessi_y_helioprojective_to_pixel(header_center, value, multiplier, additional_value):
    return _r_hp_to_px_y(header_center, value, multiplier, additional_value)


def x_helioprojective_to_pixel(value, additional_value):
    return _x_hp_to_px(value, additional_value, pixel_arcsec_x_list_placeholder, sun_center_x_list_placeholder)


def y_helioprojective_to_pixel(value, additional_value):
    return _y_hp_to_px(value, additional_value, pixel_arcsec_y_list_placeholder, sun_center_y_list_placeholder)


def x_pixel_to_helioprojective(value):
    return _x_px_to_hp(value, sun_center_x_list_placeholder, pixel_arcsec_x_list_placeholder)


def y_pixel_to_helioprojective(value):
    return _y_px_to_hp(value, sun_center_y_list_placeholder, pixel_arcsec_y_list_placeholder)


def helioprojective_to_heliographic(x, y, time):
    return _hp_to_hg(x, y, time)


def helioprojective_to_heliographic(x, y, time):
    return _hp_to_hg(x, y, time)


def find_files(url):
    soup = BeautifulSoup(requests.get(url).text, features="lxml")

    hrefs = []

    for a in soup.find_all('a'):
        hrefs.append(a['href'])

    return hrefs


def get_image_urls(dates_list):
    return _get_image_urls(dates_list)


def create_date_list(sdate, edate):
    return _create_date_list(sdate, edate)


def to_datetime_object(date_string, date_format):
    return _to_datetime_object(date_string, date_format)


def convert_timedelta(duration):
    return _convert_timedelta(duration)


import numpy as np


def mask_saturation(rast, nodata=-9999):
    return _mask_saturation(rast, nodata)


def selected_RHESSI_keys(flare_list, selected_start_date, selected_end_date):
    left_index = bisect.bisect(flare_list, selected_start_date) - 1

    right_index = bisect.bisect(flare_list, selected_end_date)

    selected_solar_flares_keys_list = flare_list[left_index:right_index]

    return selected_solar_flares_keys_list


def findpeaks(series, DELTA):
    """
    Finds extrema in a pandas series data.

    Parameters
    ----------
    series : `pandas.Series`
        The data series from which we need to find extrema.

    DELTA : `float`
        The minimum difference between data values that defines a peak.

    Returns
    -------
    minpeaks, maxpeaks : `list`
        Lists consisting of pos, val pairs for both local minima points and
        local maxima points.
    """
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
    print("maxpeaks")
    print(maxpeaks)
    return minpeaks, maxpeaks


class MathTextSciFormatter(mticker.Formatter):
    def __init__(self, fmt="%1.2e"):
        self.fmt = fmt

    def __call__(self, x, pos=None):
        s = self.fmt % x
        decimal_point = '.'
        positive_sign = '+'
        tup = s.split('e')
        significand = tup[0].rstrip(decimal_point)
        sign = tup[1][0].replace(positive_sign, '')
        exponent = tup[1][1:].lstrip('0')
        if exponent:
            exponent = '10^{%s%s}' % (sign, exponent)
        if significand and exponent:
            s = r'%s{\times}%s' % (significand, exponent)
        else:
            s = r'%s%s' % (significand, exponent)
        return "${}$".format(s)


class draggableline:
    def __init__(self, ax, XorY):
        self.ax = ax
        self.c = ax.get_figure().canvas

        self.XorY = XorY

        x = [XorY, XorY]
        y = [-1, 10]
        self.line = lines.Line2D(x, y, color='black', picker=5)
        self.ax.add_line(self.line)
        self.c.draw_idle()
        self.sid = self.c.mpl_connect('pick_event', self.clickonline)

    def clickonline(self, event):
        if event.artist == self.line:
            self.follower = self.c.mpl_connect("motion_notify_event", self.followmouse)
            self.releaser = self.c.mpl_connect("button_press_event", self.releaseonclick)

    def followmouse(self, event):
        self.line.set_xdata([event.xdata, event.xdata])
        self.c.draw_idle()

    def releaseonclick(self, event):
        self.XorY = self.line.get_xdata()[0]

        self.c.mpl_disconnect(self.releaser)
        self.c.mpl_disconnect(self.follower)


class Player(FuncAnimation):
    def __init__(self, fig, func, frames=None, init_func=None, fargs=None,
                 save_count=None, mini=0, maxi=100, pos=(0.125, 0.92), **kwargs):
        self.i = 0
        self.min = mini
        self.max = maxi
        self.runs = False
        self.forwards = True
        self.fig = fig
        self.func = func
        FuncAnimation.__init__(self, self.fig, self.func, frames=self.play(),
                               init_func=init_func, fargs=fargs,
                               save_count=save_count, **kwargs)

    def play(self):
        while self.runs:
            self.i = self.i + self.forwards - (not self.forwards)
            if self.i > self.min and self.i <= self.max:
                yield self.i

            else:
                print(self.i)
                print("return to min")
                self.i = self.min
                yield self.i

    def reset(self):
        self.i = 0

    def start(self):
        print("Start")
        self.runs = True
        self.event_source.start()

    def stop(self, event=None):
        print("stop")
        self.runs = False
        self.event_source.stop()

    def forward(self, event=None):
        print("foward")
        self.forwards = True
        self.start()

    def backward(self, event=None):
        if self.i == 0:
            self.i = self.max
            self.forwards = False
            self.start()
        else:
            self.forwards = False
            self.start()

    def oneforward(self, event=None):
        print("one foward")
        self.forwards = True
        if self.i > self.min and self.i < self.max:
            self.i = self.i + self.forwards - (not self.forwards)
        elif self.i == self.min and self.forwards:
            self.i += 1
        elif self.i == self.max and not self.forwards:
            self.i -= 1
        else:
            self.i = self.min
        self.func(self.i)
        self.fig.canvas.draw_idle()

    def onebackward(self, event=None):
        print("one backward")
        self.forwards = False
        if self.i == 0:
            self.i = self.max
        else:
            if self.i > self.min and self.i < self.max:
                self.i = self.i + self.forwards - (not self.forwards)
            elif self.i == self.min and self.forwards:
                self.i += 1
            elif self.i == self.max and not self.forwards:
                self.i -= 1
            else:
                self.i = self.max
        self.func(self.i)
        self.fig.canvas.draw_idle()
        self.fowards = True

    def onestep(self):
        if self.i > self.min and self.i < self.max:
            self.i = self.i + self.forwards - (not self.forwards)
        elif self.i == self.min and self.forwards:
            self.i += 1
        elif self.i == self.max and not self.forwards:
            self.i -= 1
        self.func(self.i)
        self.fig.canvas.draw_idle()


def create_flux_curve(time_table, total_intensity):
    tbl_meta = {'t_key': 't_value'}
    table = Table([time_table, total_intensity], names=['time', 'intensity'], meta=tbl_meta)
    table.add_index('time')

    ts_table = ts.TimeSeries(table)

    return ts_table


def split_frames(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


def flare_plot_norm_selection(normalization, ax, image_list, colormap, extent_list):
    if normalization == "asinh":
        image = ax.imshow(image_list[0], cmap=colormap,
                          extent=extent_list, norm=ImageNormalize(stretch=AsinhStretch(
                0.01)))  # , vmin=min_value, vmax=max_value_vmax, norm=ImageNormalize(stretch=AsinhStretch(0.01)))

    else:
        image = ax.imshow(image_list[0], cmap=colormap,
                          extent=extent_list)  # , vmin=min_value, vmax=max_value_vmax)

    return image


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


not_main_finalized_frame_list = []


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


def get_data(event_beginning_time, event_peak_time, event_end_time):
    print("event beginning time")
    print(event_beginning_time, event_peak_time, event_end_time)
    date_format = "%Y-%m-%dT%H:%M:%S"
    beginning_datetime_object = to_datetime_object(str(event_beginning_time)[:-4], date_format)
    peak_datetime_object = to_datetime_object(str(event_peak_time)[:-10], date_format)
    end_datetime_object = to_datetime_object(str(event_end_time)[:-4], date_format)
    beginning_to_end_flare_time_difference = end_datetime_object - beginning_datetime_object
    beginning_to_peak_flare_time_difference = peak_datetime_object - beginning_datetime_object
    peak_to_end_flare_time_difference = end_datetime_object - peak_datetime_object

    special_time_difference = beginning_to_peak_flare_time_difference * 2 + peak_to_end_flare_time_difference * 2
    print(special_time_difference)
    event_peak_time_string = event_peak_time.replace("T", " ")
    import pickle
    try:
        rhessi_image_data = master_dict[event_peak_time_string]["image_data"]
        rhessi_image_data = np.flipud(rhessi_image_data)
        with open('samplerhessidata.pickle', 'wb') as f:
            pickle.dump(rhessi_image_data, f)

        print("rhessi complete")
        rhessi_header_x = master_dict[event_peak_time_string]["rhessi_header_x"]
        rhessi_header_y = master_dict[event_peak_time_string]["rhessi_header_y"]
        multiplier = master_dict[event_peak_time_string]["rhessi_conversion"]
        rhessi_pixel_pos = np.argwhere(rhessi_image_data == rhessi_image_data.max()) * u.pixel
        rhessi_point = rhessi_pixel_pos[int(len(rhessi_pixel_pos) / 2)]
        rhessi_center_x = int(float(str(rhessi_point[1])[:-4]))
        rhessi_center_y = int(float(str(rhessi_point[0])[:-4]))
        print("rhessi info")
        print(rhessi_point)
        print(rhessi_center_x, rhessi_center_y)

        no_image_state = False
    except KeyError:
        rhessi_image_data = np.zeros((128, 128), dtype="uint8")
        rhessi_header_x = 64.5
        rhessi_header_y = 64.5
        multiplier = 20
        rhessi_center_x = 64
        rhessi_center_y = 64
        no_image_state = True

    print("rhessi image_data")
    print(rhessi_image_data)
    print(rhessi_image_data.shape)

    cropped_rhessi_data = rhessi_image_data[rhessi_center_y - cropped_dimension:rhessi_center_y + cropped_dimension,
                          rhessi_center_x - cropped_dimension:rhessi_center_x + cropped_dimension]

    def convert_timedelta(duration):
        seconds = duration.seconds
        minutes = math.ceil(seconds / 60) + 1
        return minutes

    beginning_to_end_minutes = str(convert_timedelta(beginning_to_end_flare_time_difference))
    special_beginning_to_end_minutes = str(convert_timedelta(special_time_difference))

    beginning_time_str = parse_time(event_beginning_time)
    special_beginning_time_str = parse_time(event_beginning_time) - (beginning_to_peak_flare_time_difference * 2)
    peak_time_str = str(event_peak_time) + "/30s@12s"
    beginning_to_end_query_str = str(beginning_time_str)[:-7] + "Z/{}m@12s".format(beginning_to_end_minutes)
    print("query string")
    print(beginning_to_end_query_str)

    special_beginning_to_end_query_str = str(special_beginning_time_str)[:-7] + "Z/{}m@12s".format(
        special_beginning_to_end_minutes)

    # Choose wavelnth
    wavelnth = int(initial_wavelength.get())
    print("peak time str")
    print(peak_time_str)
    peak_query, peak_s = jsoc.query('{}[{}][WAVELNTH = {}]'.format(ds[0], peak_time_str, wavelnth),
                                    key='DATE, T_OBS, FSN, DATAMEAN, TEMPCCD, EXPTIME, CPRIX1, CPRIX2, CRVAL1, CRVAL2, CDELT1, CDELT2, CROTA2, X0, Y0',
                                    seg=['image', 'spikes'])

    query, s = jsoc.query('{}[{}][WAVELNTH = {}]'.format(ds[0], special_beginning_to_end_query_str, wavelnth),
                          key='DATE, T_OBS, FSN, DATAMEAN, TEMPCCD, EXPTIME, CRPIX1, CRPIX2, CRVAL1, CRVAL2, CDELT1, CDELT2, CROTA2, X0, Y0',
                          seg=['image', 'spikes'])

    date_list = []
    query_exposure_list = []
    degrees_rotation_list = []
    cropped_map_intensity_list = []

    sun_center_x_list_placeholder.clear()
    sun_center_y_list_placeholder.clear()
    pixel_arcsec_x_list_placeholder.clear()
    pixel_arcsec_y_list_placeholder.clear()

    sun_center_x_list_placeholder.append(query.CRPIX1[0])
    sun_center_y_list_placeholder.append(query.CRPIX2[0])
    pixel_arcsec_x_list_placeholder.append(query.CDELT1[0])
    pixel_arcsec_y_list_placeholder.append(query.CDELT2[0])

    print("datasss")
    print(sun_center_y_list_placeholder[0])
    print(pixel_arcsec_y_list_placeholder[0])
    print("CPRIX1")
    print(query.X0)
    print(query.Y0)
    for i in range(0, len(query)):
        exp_time = query.EXPTIME[i]
        date = query.T_OBS[i]
        query_exposure_list.append(exp_time)
        date_list.append(date)
        rotation = query.CROTA2[i]
        degrees_rotation_list.append(rotation)
    print("degrees rotation list")
    print(degrees_rotation_list)
    query_number_of_elements = len(query_exposure_list)

    date_for_checking = peak_query.T_OBS[0]
    finalized_frame_list.clear()
    number_of_frames = int(number_of_frames_box.get())
    if auto_verification:
        print("if auto verification")
        peak_frame_index = date_list.index(date_for_checking)
        if query_exposure_list[peak_frame_index] < 1.0:
            peak_frame_index += 1
        front_temp_frame_list = list(split_frames(list(range(0, peak_frame_index)), int(number_of_frames / 2)))
        bot_temp_frame_list = list(
            split_frames(list(range(peak_frame_index + 1, query_number_of_elements)), int(number_of_frames / 2)))
        finalized_frame_list.clear()
        frame_sort_func(front_temp_frame_list, query_exposure_list, 0)
        frame_sort_func(bot_temp_frame_list, query_exposure_list, peak_frame_index)


    else:
        finalized_frame_list.clear()
        try:
            number_of_frames = int(number_of_frames_box.get())
            print("number of frames")
            print(number_of_frames)
        except ValueError:
            print("value error")
            number_of_frames = 25
        if number_of_frames > query_number_of_elements:
            number_of_frames = query_number_of_elements

        temp_frame_list = list(split_frames(list(range(0, query_number_of_elements)), number_of_frames))
        frame_sort_func(temp_frame_list, query_exposure_list, 0)
    print("finalized frame list")
    print(finalized_frame_list)

    all_image_list_dict = {}
    cropped_image_list_dict = []
    time_list = []

    brightest_pixel_dict = {}
    exposure_times_dict = {}
    selected_time_dict = {}

    peak_image_index = int(len(query) / 2)
    dummy, peak_image = fits.open('http://jsoc.stanford.edu' + s['image'][peak_image_index])
    peak_image.verify('fix')
    peak_image_data = peak_image.data
    degrees_of_rotation = int(degrees_rotation_list[peak_image_index])
    if degrees_of_rotation < 170:
        peak_image_data = np.flipud(peak_image_data)

    # peak_image_data = scipy.ndimage.rotate(peak_image_data,degrees_of_rotation, reshape = False)

    peak_pixel_pos = np.argwhere(peak_image_data == peak_image_data.max()) * u.pixel
    peak_point = peak_pixel_pos[int(len(peak_pixel_pos) / 2)]
    peak_x = int(float(str(peak_point[0])[:-4]))
    peak_y = int(float(str(peak_point[1])[:-4]))

    frame_counter = 0
    normalization = normalization_string.get()

    download_time_list = []
    non_main_time_list = []

    for i in range(0, query_number_of_elements):
        start = lmao.time()
        if i in finalized_frame_list:
            dummy, image = fits.open('http://jsoc.stanford.edu' + s['image'][i])
            image.verify('fix')
            image_data = image.data
            image_data = np.flipud(image_data)
            time1 = query.T_OBS[i][:-4]
            exposure_time = query.EXPTIME[i]
            real_time1 = to_datetime_object(time1, date_format)

            exposure_times_dict[i] = exposure_time
            pixel_pos = np.argwhere(image_data == image_data.max()) * u.pixel
            point = pixel_pos[int(len(pixel_pos) / 2)]
            x = int(float(str(point[0])[:-4]))
            y = int(float(str(point[1])[:-4]))
            cropped_image_data = image_data[peak_x - 300:peak_x + 300, peak_y - 300:peak_y + 300]
            intensity_value = np.divide(cropped_image_data.sum(), exposure_time)
            cropped_map_intensity_list.append(intensity_value)
            time_list.append(real_time1)

            all_image_list_dict[i] = image_data / exposure_time
            image_data_normalized = image_data / exposure_time
            """
            global max_value_vmax
            if image_data_normalized.max() > max_value_vmax:
                max_value_vmax = image_data_normalized.max()
                print("max value vmax")
                print(image_data_normalized.max())
                print(max_value_vmax)
            print("exposure time")
            print(exposure_time)
            """
            low_resolution_image_data = strided_rescale(image_data_normalized, 8)
            full_image_list_dict.append(low_resolution_image_data)
            new_cropped_data = image_data_normalized[peak_x - 300:peak_x + 300, peak_y - 300:peak_y + 300]
            cropped_image_list_dict.append(new_cropped_data)
            brightest_pixel_dict[frame_counter] = [x, y]
            selected_time_dict[frame_counter] = time1
            frame_counter += 1
            end = lmao.time()
            download_time = (end - start)
            download_time_list.append(download_time)
            estimated_time = (sum(download_time_list) / len(download_time_list)) * (
                    number_of_frames - (finalized_frame_list.index(i) + 1))
            print(estimated_time)
            estimated_time_minutes = int(estimated_time / 60)
            estimated_time_seconds = estimated_time - (estimated_time_minutes * 60)
            print(
                "Finished Downloading Image Array {}, Took {} Seconds, Estimated time remaining {} Minutes And {} Seconds".format(
                    str(i), str(download_time), estimated_time_minutes, estimated_time_seconds))

    all_wavelength_dict = {}
    all_wavelength_dict.clear()

    print(finalized_frame_list)
    for wavelength in wavelength_list:
        intensity_list = []
        image_data_list = []
        intensity_list.clear()
        image_data_list.clear()
        all_wavelength_dict[wavelength] = {}

        if wavelength == wavelnth:
            pass
        else:
            query, s = jsoc.query('{}[{}][WAVELNTH = {}]'.format(ds[0], beginning_to_end_query_str, wavelength),
                                  key='DATE, T_OBS, FSN, DATAMEAN, TEMPCCD, EXPTIME, CRPIX1, CRPIX2, CRVAL1, CRVAL2, CDELT1, CDELT2, CROTA2, X0, Y0',
                                  seg=['image', 'spikes'])
            query_number_of_elements = len(query)

            finalized_frame_list.clear()
            number_of_frames = int(number_of_frames_box.get())
            if auto_verification:
                print("if auto verification")
                peak_frame_index = date_list.index(date_for_checking)
                if query_exposure_list[peak_frame_index] < 1.0:
                    peak_frame_index += 1
                front_temp_frame_list = list(split_frames(list(range(0, peak_frame_index)), int(number_of_frames / 2)))
                bot_temp_frame_list = list(
                    split_frames(list(range(peak_frame_index + 1, query_number_of_elements)),
                                 int(number_of_frames / 2)))
                finalized_frame_list.clear()
                frame_sort_func(front_temp_frame_list, query_exposure_list, 0)
                frame_sort_func(bot_temp_frame_list, query_exposure_list, peak_frame_index)
                print("finalized frame")
                print(finalized_frame_list)
            else:
                finalized_frame_list.clear()
                try:
                    number_of_frames = int(number_of_frames_box.get())
                    print("number of frames")
                    print(number_of_frames)
                except ValueError:
                    print("value error")
                    number_of_frames = 25
                if number_of_frames > query_number_of_elements:
                    number_of_frames = query_number_of_elements

                temp_frame_list = list(split_frames(list(range(0, query_number_of_elements)), number_of_frames))
                frame_sort_func(temp_frame_list, query_exposure_list, 0)
                print("finalized frame")
                print(finalized_frame_list)

            if wavelength == 335:
                for number in finalized_frame_list:
                    print("qiuerryy")
                    print(query)
                    print(finalized_frame_list)
                    time1 = query.T_OBS[number][:-4]
                    real_time1 = to_datetime_object(time1, date_format)
                    non_main_time_list.append(real_time1)

            for i in range(0, len(query)):
                if i in finalized_frame_list:
                    start = lmao.time()
                    dummy, image = fits.open('http://jsoc.stanford.edu' + s['image'][i])
                    image.verify('fix')
                    image_data = image.data
                    image_data = np.flipud(image_data)
                    exposure_time = query.EXPTIME[i]
                    cropped_image_data = image_data[peak_x - 300:peak_x + 300, peak_y - 300:peak_y + 300]
                    intensity_value = np.divide(cropped_image_data.sum(), exposure_time)
                    intensity_list.append(intensity_value)
                    image_data_list.append(image_data)

                    end = lmao.time()
                    download_time = (end - start)
                    download_time_list.append(download_time)
                    estimated_time = (sum(download_time_list) / len(download_time_list)) * (
                            number_of_frames - (finalized_frame_list.index(i) + 1))
                    print(estimated_time)
                    estimated_time_minutes = int(estimated_time / 60)
                    estimated_time_seconds = estimated_time - (estimated_time_minutes * 60)
                    print(
                        "Finished Downloading Image Array {}, Took {} Seconds, Estimated time remaining {} Minutes And {} Seconds".format(
                            str(i), str(download_time), estimated_time_minutes, estimated_time_seconds))

            all_wavelength_dict[wavelength]["intensity_list"] = intensity_list
            all_wavelength_dict[wavelength]["image_data"] = image_data_list

    print("wavelength dict")
    print(all_wavelength_dict)

    return full_image_list_dict, time_list, brightest_pixel_dict, cropped_image_list_dict, peak_x, peak_y, exposure_times_dict, all_image_list_dict, selected_time_dict, rhessi_image_data, rhessi_center_x, rhessi_center_y, cropped_rhessi_data, rhessi_header_x, rhessi_header_y, multiplier, no_image_state, cropped_map_intensity_list, all_wavelength_dict, non_main_time_list


def scale_image(input_image_path,
                output_image_path,
                width=None,
                height=None
                ):
    original_image = Image.open(input_image_path)
    w, h = original_image.size
    print('The original image size is {wide} wide x {height} '
          'high'.format(wide=w, height=h))
    if width and height:
        max_size = (width, height)
    elif width:
        max_size = (width, h)
    elif height:
        max_size = (w, height)
    else:
        # No width or height specified
        raise RuntimeError('Width or height required!')
    original_image.thumbnail(max_size, Image.ANTIALIAS)
    original_image.save(output_image_path)
    scaled_image = Image.open(output_image_path)
    width, height = scaled_image.size
    print('The scaled image size is {wide} wide x {height} '
          'high'.format(wide=width, height=height))


image1 = requests.get("https://i.imgur.com/klMBNo3.jpg")

file = open("lockheedlogo1.jpg", "wb")
file.write(image1.content)
file.close()

image2 = requests.get("https://i.imgur.com/pY1zShM.png")

file = open("lockheedlogo2.png", "wb")
file.write(image2.content)
file.close()

path = "lockheedfinal1.jpg"
path2 = "lockheedfinal2.png"

scale_image(input_image_path='lockheedlogo1.jpg',
            output_image_path='lockheedfinal1.jpg',
            width=100)

scale_image(input_image_path='lockheedlogo2.png',
            output_image_path='lockheedfinal2.png',
            width=100)

# Creates a Tkinter-compatible photo image, which can be used everywhere Tkinter expects an image object.
img = ImageTk.PhotoImage(Image.open(path))
img3 = ImageTk.PhotoImage(Image.open(path2))

# The Label widget is a standard Tkinter widget used to display a text or image on the screen.
panel = Label(root, bg="white", image=img)
panel1 = Label(root, bg="white", image=img3)

# The Pack geometry manager packs widgets in rows or columns.
panel.place(x=80, y=230)
panel1.place(x=220, y=230)

start_date = ''
end_date = ''

start_date_entry_label = Label(root, bg="white", text="Start Date:(yyyy-mm-dd hr-min) ")
end_date_entry_label = Label(root, bg="white", text="End Date:(yyyy-mm-dd hr-min) ")

start_date_entry_label.place(x=20, y=30)
end_date_entry_label.place(x=20, y=70)

end_date_string = StringVar(root)
start_date_entry_box = Entry(root, width=14, highlightbackground="white")
end_date_entry_box = Entry(root, textvariable=end_date_string, width=14, highlightbackground="white")

start_date_entry_box.place(x=250, y=27)
end_date_entry_box.place(x=250, y=67)

goes_minimum_value_entry_label = Label(root, bg="white", text="Minimum Flare Class (Ex. M2.2)")
goes_maximum_value_entry_label = Label(root, bg="white", text="Maximum Flare Class (Ex. M2.2)")

goes_minimum_value_entry_label.place(x=20, y=110)
goes_maximum_value_entry_label.place(x=20, y=150)

goes_minimum_value_entry_box = Entry(root, highlightbackground="white", width=14)
goes_maximum_value_entry_box = Entry(root, highlightbackground="white", width=14)

goes_minimum_value_entry_box.insert(0, 'C1')
goes_maximum_value_entry_box.insert(0, 'None')

goes_minimum_value_entry_box.place(x=250, y=107)
goes_maximum_value_entry_box.place(x=250, y=147)

past_input_dict = {}

cmap_list = ["sdoaia94", "sdoaia131", "sdoaia171", "sdoaia193", "sdoaia211", "sdoaia304", "sdoaia335"]


def flares_event_dropdown_function():
    global flare_selected
    global one_frame_state
    global graph_type_regular_bool
    global graph_type_changed_bool
    global GOES_frame
    global slider_changed
    global user_clicked_map_bool
    global flare_selected
    flare_selected = False
    one_frame_state = False
    graph_type_regular_bool = True
    graph_type_changed_bool = False
    GOES_frame = True
    slider_changed = False
    animation_interval_changed = False
    user_clicked_map_bool = False
    slider_used_func = False
    flare_selected = True
    print("getting flare")
    # canvas11.get_tk_widget().destroy()
    # frame2.destroy()
    # create_frame2()
    change_frame_func()

    change_frame_button.place(x=690 + space, y=652)

    flares_event_dropdown_button.place(x=810 + space, y=380)

    flare_event_name = initial_menu_value.get()
    index = flares_dropdown_menu.index(flare_event_name)

    event_beginning_time1 = str(manual_selected_time_range[0]).replace(" ", 'T')
    event_end_time1 = str(manual_selected_time_range[1]).replace(" ", 'T')
    event_peak_time1 = str(real_peak_time_list[index]).replace(" ", 'T')

    print("dates")
    print(event_beginning_time1)

    wavelnth = int(initial_wavelength.get())
    wavelength_index = wavelength_list.index(wavelnth)
    colormap = cmap_list[wavelength_index]

    full_image_dict, time_list, brightest_pixel_dict, cropped_image_list, x_center, y_center, exposure_dict, all_images_dict, chosen_time_dict, rhessi_image_data, rhessi_x_center, rhessi_y_center, cropped_rhessi_data, rhessi_header_x, rhessi_header_y, rhessi_multiplier, no_image_state, cropped_intensity_list, all_wavelength_dict, non_main_time_list = get_data(
        event_beginning_time1, event_peak_time1, event_end_time1)
    print("timessss")
    print(event_beginning_time1, event_peak_time1, event_end_time1)

    """
    json_list = []
    for data in full_image_dict:
        print("here")
        json_list.append(data.tolist())

    json.dump(json_list, codecs.open("fulldict6.json", 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4) ### this saves the array in .json format
    """

    fig1, ax16 = plt.subplots(figsize=(3.44, 4))
    fig2, ax12 = plt.subplots(figsize=(3.44, 3.44))
    fig3, ax13 = plt.subplots(figsize=(3.42, 3.42))
    fig4, ax14 = plt.subplots(figsize=(5, 4.01))
    fig5, ax15 = plt.subplots(figsize=(3.33, 4))
    fig6, ax30 = plt.subplots(figsize=(3.33, 4))
    # fig7, ax30 = plt.subplots(figsize=(0.1,0.1))
    # fig8, ax18 = plt.subplots(figsize=(2, 3.42))

    ax15.set_visible(False)
    ax30.set_visible(False)
    # ax18.set_visible(False)

    canvas12 = FigureCanvasTkAgg(fig2, frame2)
    canvas13 = FigureCanvasTkAgg(fig3, frame2)
    canvas14 = FigureCanvasTkAgg(fig4, frame2)
    canvas16 = FigureCanvasTkAgg(fig1, frame2)
    canvas17 = FigureCanvasTkAgg(fig5, frame2)
    canvas18 = FigureCanvasTkAgg(fig6, frame2)

    # canvas20 = FigureCanvasTkAgg(fig8, frame2)

    canvas12.get_tk_widget().grid(row=1, column=0)
    canvas14.get_tk_widget().grid(row=0, column=1)
    canvas13.get_tk_widget().grid(row=1, column=1)
    canvas16.get_tk_widget().grid(row=0, column=0)
    canvas18.get_tk_widget().grid(row=2, column=1)
    canvas17.get_tk_widget().grid(row=2, column=0)
    # canvas20.get_tk_widget().grid(row=1, column=2)

    canvas13.get_tk_widget().place(x=340, y=401)

    try:
        minimum_value = str(goes_minimum_value_entry_box.get())
    except ValueError:
        minimum_value = ""
    try:
        maximum_value = str(goes_maximum_value_entry_box.get())
    except ValueError:
        maximum_value = ""

    if minimum_value == "None":
        minimum_value = ""
    if maximum_value == "None":
        maximum_value = ""

    flares_for_rhessi_comparison.clear()

    for i in range(0, len(flare_class_list)):
        event_peak_time3 = peak_time_list[i]
        event_beginning_time3 = start_time_list[i]

        flare_class = flare_class_list[i]
        # ax11.legend(loc=2)

        beginning_time_data = str(event_beginning_time3)[:-7]
        time_in_string = str(event_peak_time3)[:19]

        maximum_point = max(goes.data['xrsb'][time_in_string])
        maximum_point_list.append(maximum_point)

    hover_status = ["yes"]

    number = [initial_index]

    def normalized_plot(plot_list):
        new_list = []
        max_value = max(plot_list)
        for intensity in plot_list:
            new_val = intensity / max_value
            new_list.append(new_val)
        return new_list

    square_side_length = 300
    flare_text = StringVar()
    user_flare_text = StringVar()

    global flare_center_coordinates

    user_flare_coordinates = Label(root, bg="white", textvariable=user_flare_text)
    user_flare_coordinates.place(x=820 + space, y=625)
    flare_center_coordinates = Label(root, bg="white", textvariable=flare_text)
    flare_center_coordinates.place(x=820 + space, y=580)

    coordinate_system_string = StringVar()
    coordinate_system_string.set("helioprojective")
    coordinate_system_list_placeholder = ["helioprojective"]

    def change_coordinate_system():
        global coord_type_changed_bool
        coord_type_changed_bool = True
        coordinate_system_list_placeholder.clear()
        coordinate_system_list_placeholder.append(coordinate_system_string.get())

    coord_button_placement_y = 530
    helioprojective_radio_button = Radiobutton(root, highlightbackground="white", bg="white", text='Helioprojective',
                                               variable=coordinate_system_string,
                                               value='helioprojective', command=change_coordinate_system)
    helioprojective_radio_button.place(x=690 + space, y=coord_button_placement_y + 40)
    heliographical_radio_button = Radiobutton(root, highlightbackground="white", bg="white", text='Heliographical',
                                              variable=coordinate_system_string,
                                              value='heliographical', command=change_coordinate_system)
    heliographical_radio_button.place(x=690 + space, y=coord_button_placement_y + 62)
    pixel_radio_button = Radiobutton(root, highlightbackground="white", text='Pixel', bg="white",
                                     variable=coordinate_system_string, value='pixel',
                                     command=change_coordinate_system)
    pixel_radio_button.place(x=690 + space, y=coord_button_placement_y + 84)

    center_display_method_label = Label(root, bg="white", text="Coordinate\nSystem", justify=CENTER)
    center_display_method_label.place(x=710 + space, y=coord_button_placement_y)

    random_time = time_list[0]

    rhessi_helioprojective_x = rhessi_pixel_to_helioprojective_x(rhessi_header_x, rhessi_x_center, rhessi_multiplier)
    rhessi_helioprojective_y = -rhessi_pixel_to_helioprojective_x(rhessi_header_y, rhessi_y_center, rhessi_multiplier)
    rhessi_heliographical_coords = helioprojective_to_heliographic(rhessi_helioprojective_x, rhessi_helioprojective_y,
                                                                   random_time)

    rhessi_axis_dict.clear()

    rhessi_axis_dict["rhessi_helioprojective_x"] = rhessi_helioprojective_x
    rhessi_axis_dict["rhessi_helioprojective_y"] = rhessi_helioprojective_y
    rhessi_axis_dict["rhessi_heliographical_coords"] = rhessi_heliographical_coords
    rhessi_axis_dict["rhessi_center_x"] = rhessi_x_center
    rhessi_axis_dict["rhessi_center_y"] = rhessi_y_center

    first_rhessi_x_center = rhessi_x_center
    first_rhessi_y_center = rhessi_y_center
    first_rhessi_helioprojective_x = rhessi_helioprojective_x
    first_rhessi_helioprojective_y = rhessi_helioprojective_y
    first_rhessi_heliographical_coords = rhessi_heliographical_coords

    rhessi_flare_text = StringVar()
    user_rhessi_flare_text = StringVar()

    def coord_change_func():
        coord_type = coordinate_system_string.get()
        rhessi_flare_text.set("")
        user_rhessi_flare_text.set("")
        if coord_type == "helioprojective":
            rhessi_flare_text.set("({},{})".format(first_rhessi_helioprojective_x, first_rhessi_helioprojective_y
                                                   ))
        if coord_type == "heliographical":
            rhessi_flare_text.set("({})".format(first_rhessi_heliographical_coords))
        if coord_type == "pixel":
            rhessi_flare_text.set("({},{})".format(int(first_rhessi_x_center) * 32,
                                                   int(first_rhessi_y_center) * 32))

        if coord_type == "helioprojective":
            user_rhessi_flare_text.set("({},{})".format(rhessi_axis_dict["rhessi_helioprojective_x"],
                                                        rhessi_axis_dict["rhessi_helioprojective_y"]))
        if coord_type == "heliographical":
            user_rhessi_flare_text.set("({})".format(rhessi_axis_dict["rhessi_heliographical_coords"]))
        if coord_type == "pixel":
            user_rhessi_flare_text.set("({},{})".format(int(rhessi_axis_dict["rhessi_center_x"]) * 32,
                                                        int(rhessi_axis_dict["rhessi_center_y"]) * 32))

    coord_change_func()

    number_of_images = len(cropped_image_list) - 1

    sun_center_x = sun_center_x_list_placeholder[0]
    sun_center_y = sun_center_y_list_placeholder[0]
    pixel_arcsec_x = pixel_arcsec_x_list_placeholder[0]
    pixel_arcsec_y = pixel_arcsec_y_list_placeholder[0]

    original_solar_movies_dict = {}
    changing_solar_movies_dict = {}
    key.append(0)

    def orignal_solar_dict_function():
        """
        original_rect = patches.Rectangle((y_center/8 - square_side_length/8, x_center/8 - square_side_length/8),
                                     square_side_length /4,
                                         square_side_length /4, linewidth=1.5, edgecolor='', facecolor='none')


        original_solar_movies_dict["rectangle"] = original_rect
        """

        for i in range(0, number_of_images + 1):
            flare_x_center = brightest_pixel_dict[i][1]
            flare_y_center = brightest_pixel_dict[i][0]
            helioprojective_center_x = x_pixel_to_helioprojective(flare_x_center)
            helioprojective_center_y = y_pixel_to_helioprojective(flare_y_center)
            heliographical_center_coords = helioprojective_to_heliographic(helioprojective_center_x,
                                                                           helioprojective_center_y, random_time)

            original_solar_movies_dict[i] = [flare_x_center, flare_y_center, helioprojective_center_x,
                                             helioprojective_center_y, heliographical_center_coords]

        # have to flip here
        small_x = y_center - square_side_length
        small_y = x_center - square_side_length
        big_x = y_center + square_side_length
        big_y = x_center + square_side_length

        helioprojective_small_x = x_pixel_to_helioprojective(small_x)
        helioprojective_small_y = y_pixel_to_helioprojective(small_y)
        rounded_small_x = roundup(helioprojective_small_x)
        rounded_small_y = rounddown(helioprojective_small_y)

        print("axis labels")
        print(small_y)
        print(helioprojective_small_y)
        print(rounded_small_y)
        first_x = x_helioprojective_to_pixel(rounded_small_x, 0)
        second_x = x_helioprojective_to_pixel(rounded_small_x, 100)
        third_x = x_helioprojective_to_pixel(rounded_small_x, 200)
        first_y = y_helioprojective_to_pixel(rounded_small_y, 0)
        second_y = y_helioprojective_to_pixel(rounded_small_y, 100)
        third_y = y_helioprojective_to_pixel(rounded_small_y, 200)

        changing_solar_movies_dict["helioprojective_small_x"] = helioprojective_small_x
        changing_solar_movies_dict["helioprojective_small_y"] = helioprojective_small_y
        changing_solar_movies_dict["small_x"] = small_x
        changing_solar_movies_dict["small_y"] = small_y
        changing_solar_movies_dict["big_x"] = big_x
        changing_solar_movies_dict["big_y"] = big_y
        changing_solar_movies_dict["rounded_small_x"] = rounded_small_x
        changing_solar_movies_dict["rounded_small_y"] = rounded_small_y
        changing_solar_movies_dict["first_x"] = first_x
        changing_solar_movies_dict["second_x"] = second_x
        changing_solar_movies_dict["third_x"] = third_x
        changing_solar_movies_dict["first_y"] = first_y
        changing_solar_movies_dict["second_y"] = second_y
        changing_solar_movies_dict["third_y"] = third_y
        changing_solar_movies_dict["user_selected_center_pixel"] = y_center, x_center
        user_projective_x = x_pixel_to_helioprojective(y_center)
        user_projective_y = y_pixel_to_helioprojective(x_center)
        changing_solar_movies_dict["user_selected_center_helioprojective"] = user_projective_x, user_projective_y
        changing_solar_movies_dict["user_selected_center_heliographical"] = helioprojective_to_heliographic(
            user_projective_x, user_projective_y, random_time)

    orignal_solar_dict_function()

    # ax12.set_xticks([sun_center_x-(1000/pixel_arcsec_x), sun_center_x-(500/pixel_arcsec_x), sun_center_x, sun_center_x+(500/pixel_arcsec_x), sun_center_x+(1000/pixel_arcsec_x)])

    small_x = changing_solar_movies_dict["small_x"]
    small_y = changing_solar_movies_dict["small_y"]
    big_x = changing_solar_movies_dict["big_x"]
    big_y = changing_solar_movies_dict["big_y"]

    norm = normalization_string.get()
    global im, im2
    im = flare_plot_norm_selection(norm, ax12, full_image_dict, colormap, [0, 512, 512, 0])
    im2 = flare_plot_norm_selection(norm, ax13, cropped_image_list, colormap, [small_x, big_x, big_y, small_y])
    im3 = ax16.imshow(rhessi_image_data, cmap='rhessi', extent=[0, 128, 128, 0])
    rhessi_center_point_marker = ax16.plot(rhessi_x_center, rhessi_y_center, 'wx', color='green', fillstyle='none',
                                           markersize=15)

    initial_x = 12.5
    adding_x = 25
    ax16.set_xticks(
        [initial_x, initial_x + adding_x, initial_x + 2 * adding_x, initial_x + 3 * adding_x, initial_x + 4 * adding_x])
    ax16.set_xticklabels(['-1000"', '-500"', '0"', '500"', '1000"'])

    ax16.set_yticks(
        [initial_x, initial_x + adding_x, initial_x + 2 * adding_x, initial_x + 3 * adding_x, initial_x + 4 * adding_x])
    ax16.set_yticklabels(['1000"', '500"', '0"', '-500"', '-1000"'])

    flare_center_point, = ax13.plot([], [], 'wx', color='red', fillstyle='none', markersize=15)
    flare_rect, = ax12.plot([], [], 's', color='green', fillstyle='none', markersize=30)
    user_flare_rect, = ax12.plot([], [], 's', color='red', fillstyle='none', markersize=30)

    ax12.set_xticks(
        [(sun_center_x - (1000 / pixel_arcsec_x)) / 8, (sun_center_x - (500 / pixel_arcsec_x)) / 8,
         sun_center_x / 8,
         (sun_center_x + (500 / pixel_arcsec_x)) / 8, (sun_center_x + (1000 / pixel_arcsec_x)) / 8])
    ax12.set_xticklabels(['-1000"', '-500"', '0"', '500"', '1000"'])
    # ax12.set_yticks([sun_center_y-(1000/pixel_arcsec_y), sun_center_y-(500/pixel_arcsec_y), sun_center_y, sun_center_y+(500/pixel_arcsec_y), sun_center_y+(1000/pixel_arcsec_x)])
    ax12.set_yticks(
        [(sun_center_y - (1000 / pixel_arcsec_y)) / 8, (sun_center_y - (500 / pixel_arcsec_y)) / 8,
         sun_center_y / 8,
         (sun_center_y + (500 / pixel_arcsec_y)) / 8, (sun_center_y + (1000 / pixel_arcsec_x)) / 8])
    ax12.set_yticklabels(['1000"', '500"', '0"', '-500"', '-1000"'])
    # flare_x_center_, flare_y_center_, helioprojective_flare_center_x_, helioprojective_flare_center_y_, original_rect_one = \
    #   original_solar_movies_dict[0]

    big_plot_time_string = StringVar()
    big_plot_time_label = Label(frame2, bg="white", text="fdsfdsfdsfdsfdsfdsfds", justify=CENTER,
                                textvariable=big_plot_time_string, font=(None, 16))
    big_plot_time_label.place(x=100, y=400)

    def update(i):
        if key:
            im.set_array(full_image_dict[current_frame_number])
            print("time printed")

            try:
                flare_x_center, flare_y_center, helioprojective_flare_center_x, helioprojective_flare_center_y, heliographic_coords = \
                    original_solar_movies_dict[current_frame_number - 1]
            except KeyError:
                flare_x_center, flare_y_center, helioprojective_flare_center_x, helioprojective_flare_center_y, heliographic_coords = \
                    original_solar_movies_dict[0]


        else:
            im.set_array(full_image_dict[current_frame_number])

            user_flare_rect.set_data(changing_solar_movies_dict["ix"], changing_solar_movies_dict["iy"])
            flare_x_center, flare_y_center, helioprojective_flare_center_x, helioprojective_flare_center_y, heliographic_coords = \
                original_solar_movies_dict[current_frame_number]

        fig2.suptitle(" ")
        flare_text.set("")
        user_flare_text.set("")
        if coordinate_system_list_placeholder[0] == "helioprojective":
            flare_text.set(
                "({},{})".format("%.0f" % helioprojective_flare_center_x, "%.0f" % helioprojective_flare_center_y))
            user_flare_text.set(
                "({},{})".format("%.0f" % changing_solar_movies_dict["user_selected_center_helioprojective"][0],
                                 "%.0f" % (changing_solar_movies_dict["user_selected_center_helioprojective"][1])))
        elif coordinate_system_list_placeholder[0] == "heliographical":
            flare_text.set("({})".format(
                heliographic_coords))
            """"
                helioprojective_to_heliographic(helioprojective_flare_center_x, helioprojective_flare_center_y,
                                                random_time)))
            """
            print("user selected heliographical")
            """
            user_flare_text.set("({})".format(
                helioprojective_to_heliographic(changing_solar_movies_dict["user_selected_center_helioprojective"][0],
                                                changing_solar_movies_dict["user_selected_center_helioprojective"][1],
                                                random_time)))
            """
            user_flare_text.set("({})".format(
                changing_solar_movies_dict["user_selected_center_heliographical"]))
        else:
            flare_text.set("({},{})".format(flare_x_center, flare_y_center))
            user_flare_text.set("({},{})".format(changing_solar_movies_dict["user_selected_center_pixel"][0],
                                                 (changing_solar_movies_dict["user_selected_center_pixel"][1])))

        flare_rect.set_data(y_center / 8, x_center / 8)

        if key:
            return [im, flare_rect, ]
        else:
            return [im, flare_rect, user_flare_rect, ]

    def change_ax_13_axis():
        helioprojective_small_x = changing_solar_movies_dict["helioprojective_small_x"]
        helioprojective_small_y = changing_solar_movies_dict["helioprojective_small_y"]
        rounded_small_x = changing_solar_movies_dict["rounded_small_x"]
        rounded_small_y = changing_solar_movies_dict["rounded_small_y"]
        first_x = changing_solar_movies_dict["first_x"]
        second_x = changing_solar_movies_dict["second_x"]
        third_x = changing_solar_movies_dict["third_x"]
        first_y = changing_solar_movies_dict["first_y"]
        second_y = changing_solar_movies_dict["second_y"]
        third_y = changing_solar_movies_dict["third_y"]
        print("axis edges")
        print(first_y)
        print(second_y)
        print(third_y)
        print(rounded_small_y)
        print((
            [str(rounded_small_y) + '"', str(rounded_small_y - 100) + '"', str(rounded_small_y - 200) + '"',
             str(rounded_small_y - 300) + '"']))

        if (rounded_small_x - helioprojective_small_x) < 59:
            fourth_x = x_helioprojective_to_pixel(rounded_small_x, 300)
            ax13.set_xticks([first_x, second_x, third_x, fourth_x])
            ax13.set_xticklabels(
                [str(rounded_small_x) + '"', str(rounded_small_x + 100) + '"', str(rounded_small_x + 200) + '"',
                 str(rounded_small_x + 300) + '"'])
        else:
            ax13.set_xticks([first_x, second_x, third_x])
            ax13.set_xticklabels(
                [str(rounded_small_x) + '"', str(rounded_small_x + 100) + '"', str(rounded_small_x + 200) + '"'])

        if (rounded_small_y - helioprojective_small_y) < 1000:

            fourth_y = y_helioprojective_to_pixel(rounded_small_y, 300)
            ax13.set_yticks([first_y, second_y, third_y, fourth_y])
            ax13.set_yticklabels(
                [str(rounded_small_y) + '"', str(rounded_small_y - 100) + '"', str(rounded_small_y - 200) + '"',
                 str(rounded_small_y - 300) + '"'])
        else:
            ax13.set_yticks([first_y, second_y, third_y])
            ax13.set_yticklabels(
                [str(rounded_small_y) + '"', str(rounded_small_y - 100) + '"', str(rounded_small_y - 200) + '"'])

    change_ax_13_axis()

    def update2(i):
        global animation_interval_changed

        im2.set_array(cropped_image_list[current_frame_number])

        flare_x_center = original_solar_movies_dict[current_frame_number][0]
        flare_y_center = original_solar_movies_dict[current_frame_number][1]
        if (small_x < flare_x_center < big_x) and (
                small_y < flare_y_center < big_y):
            flare_center_point.set_data(flare_x_center, flare_y_center)

        if animation_interval_changed == True:
            global animation_interval
            big_plot._interval = animation_interval
            small_plot._interval = animation_interval
            graph_plot._interval = animation_interval
            animation_interval_changed = False

        # fig3.tight_layout()

        return [im2, flare_center_point, ]

    ax17 = ax14.twinx()

    canvas13.draw()
    canvas12.draw()
    fig3.tight_layout()
    fig2.tight_layout()
    fig1.tight_layout()

    graph_type_string = StringVar()
    graph_type_string.set("regular")

    def graph_type_func():
        global graph_type_regular_bool
        global graph_type_changed_bool
        if graph_type_string.get() == "regular":
            graph_type_regular_bool = True
        else:
            graph_type_regular_bool = False
        graph_type_changed_bool = True

    normalized_cropped_intensity_list = normalized_plot(cropped_intensity_list)

    user_cropped_map_intensity_list = []

    regular_button = Radiobutton(root, highlightbackground="white", bg="white", text='Regular',
                                 variable=graph_type_string, command=graph_type_func,
                                 value='regular')

    normalized_button = Radiobutton(root, highlightbackground="white", bg="white", text='Normalized',
                                    variable=graph_type_string, command=graph_type_func,
                                    value='normalized')

    euv_graph_label_placements = 445
    regular_button.place(x=690 + space, y=euv_graph_label_placements + 35)
    normalized_button.place(x=690 + space, y=euv_graph_label_placements + 57)
    euv_graph_type_label = Label(root, bg="white", text="EUV Graph\nType", justify=CENTER)
    euv_graph_type_label.place(x=705 + space, y=euv_graph_label_placements)

    user_euv_plot_list = []

    only_flare_plot = goes.truncate(TimeRange(event_beginning_time1, event_end_time1))
    only_flare_plot.plot(axes=ax14)

    print("listssss")
    print(time_list)
    print(cropped_intensity_list)
    print(all_wavelength_dict[94]['intensity_list'])

    euv_flux_plot, = ax17.plot_date(time_list, cropped_intensity_list, marker='', linestyle='-', color='blue',
                                    label="WL 171")
    print("nba")
    print(len(non_main_time_list))
    print(non_main_time_list)
    print(len(all_wavelength_dict[94]))
    print(all_wavelength_dict[94])

    euv_flux_plot_94, = ax17.plot_date(non_main_time_list, all_wavelength_dict[94]['intensity_list'], marker='',
                                       linestyle='-',
                                       color='brown',
                                       label="WL 94")

    euv_flux_plot_131, = ax17.plot_date(non_main_time_list, all_wavelength_dict[131]['intensity_list'], marker='',
                                        linestyle='-', color='black',
                                        label="WL 131")
    euv_flux_plot_193, = ax17.plot_date(non_main_time_list, all_wavelength_dict[193]['intensity_list'], marker='',
                                        linestyle='-', color='pink',
                                        label="WL 193")
    euv_flux_plot_211, = ax17.plot_date(non_main_time_list, all_wavelength_dict[211]['intensity_list'], marker='',
                                        linestyle='-', color='grey',
                                        label="WL 211")
    euv_flux_plot_304, = ax17.plot_date(non_main_time_list, all_wavelength_dict[304]['intensity_list'], marker='',
                                        linestyle='-', color='purple',
                                        label="WL 304")
    euv_flux_plot_335, = ax17.plot_date(non_main_time_list, all_wavelength_dict[335]['intensity_list'], marker='',
                                        linestyle='-', color='green',
                                        label="WL 335")

    normalized_euv_plot_94 = normalized_plot(all_wavelength_dict[94]['intensity_list'])
    normalized_euv_plot_131 = normalized_plot(all_wavelength_dict[131]['intensity_list'])
    normalized_euv_plot_193 = normalized_plot(all_wavelength_dict[193]['intensity_list'])
    normalized_euv_plot_211 = normalized_plot(all_wavelength_dict[211]['intensity_list'])
    normalized_euv_plot_304 = normalized_plot(all_wavelength_dict[304]['intensity_list'])
    normalized_euv_plot_335 = normalized_plot(all_wavelength_dict[335]['intensity_list'])

    user_euv_flux_plot, = ax17.plot_date([], [], marker='', linestyle='-', color='red', label="User Intensity")

    user_euv_plot_list.append(user_euv_flux_plot)

    h2, l2 = ax17.get_legend_handles_labels()
    h1, l1 = ax14.get_legend_handles_labels()
    ax14.legend(h1 + h2, l1 + l2, loc=2)
    ax17.legend().set_visible(False)

    ax17.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
    ax14.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
    ax14.set_xlabel('Time (dd-hh-min)')
    ax14.set_ylabel('$Watts/m^{-2}$')
    ax17.set_ylabel('EUV Flux (DN/sec)')

    fig4.autofmt_xdate()
    date_format = "%Y-%m-%dT%H:%M:%S"
    time_placeholder = chosen_time_dict[0]
    real_time_placeholder = to_datetime_object(time_placeholder, date_format)
    vertical_line = ax14.axvline(real_time_placeholder, ls='-', color='black', lw=1, zorder=10)

    ax17.relim()
    ax17.autoscale_view()
    canvas14.draw()

    def update3(i):
        global current_frame_number
        global coord_type_changed_bool
        print("one frame state")
        print(one_frame_state)
        print("not frame state")
        current_frame_number += 1
        slider_val.set(current_frame_number)

        if coord_type_changed_bool:
            coord_change_func()
            coord_type_changed_bool = False

        if current_frame_number > number_of_images or current_frame_number < 0:
            current_frame_number = 0

        if len(ax17.lines) > 9:
            ax17.lines[9].remove()
        global graph_type_regular_bool
        global graph_type_changed_bool
        global user_clicked_map_bool
        time = chosen_time_dict[current_frame_number]
        real_time = to_datetime_object(time, date_format)
        vertical_line.set_xdata(real_time)
        time_string = str(time).replace("T", " ")
        big_plot_time_string.set(time_string)

        if graph_type_changed_bool:
            if graph_type_regular_bool:
                euv_flux_plot.set_data(time_list, cropped_intensity_list)
                ax17.relim()
                euv_flux_plot_94.set_data(non_main_time_list, normalized_euv_plot_94)
                ax17.relim()
                euv_flux_plot_131.set_data(non_main_time_list, normalized_euv_plot_131)
                ax17.relim()
                euv_flux_plot_193.set_data(non_main_time_list, normalized_euv_plot_193)
                ax17.relim()
                euv_flux_plot_211.set_data(non_main_time_list, normalized_euv_plot_211)
                ax17.relim()
                euv_flux_plot_304.set_data(non_main_time_list, normalized_euv_plot_304)
                ax17.relim()
                euv_flux_plot_335.set_data(non_main_time_list, normalized_euv_plot_335)
                ax17.relim()
                user_euv_flux_plot.set_data(non_main_time_list, user_cropped_map_intensity_list)


            else:
                euv_flux_plot.set_data(time_list, normalized_cropped_intensity_list)
                ax17.relim()
                euv_flux_plot_94.set_data(non_main_time_list, all_wavelength_dict[94]['intensity_list'])
                ax17.relim()
                euv_flux_plot_131.set_data(non_main_time_list, all_wavelength_dict[131]['intensity_list'])
                ax17.relim()
                euv_flux_plot_193.set_data(non_main_time_list, all_wavelength_dict[193]['intensity_list'])
                ax17.relim()
                euv_flux_plot_211.set_data(non_main_time_list, all_wavelength_dict[211]['intensity_list'])
                ax17.relim()
                euv_flux_plot_304.set_data(non_main_time_list, all_wavelength_dict[304]['intensity_list'])
                ax17.relim()
                euv_flux_plot_335.set_data(non_main_time_list, all_wavelength_dict[335]['intensity_list'])
                ax17.relim()
                try:
                    user_euv_flux_plot.set_data(time_list, changing_solar_movies_dict[
                        "user_selected_normalized_graph"])
                except KeyError:
                    pass

            print("redrew canvas")
            ax17.relim()
            ax17.autoscale_view()
            canvas14.draw()
            graph_type_changed_bool = False  # test, change if nessacary

        # fig4.tight_layout()

        if user_clicked_map_bool:
            if graph_type_regular_bool:
                user_euv_flux_plot.set_data(time_list, user_cropped_map_intensity_list)
            else:
                user_euv_flux_plot.set_data(time_list, changing_solar_movies_dict[
                    "user_selected_normalized_graph"])
            print("redrew canvas2")
            ax17.relim()
            ax17.autoscale_view()
            canvas14.draw()
            user_clicked_map_bool = False
        print("comes here")
        return vertical_line, user_euv_flux_plot, euv_flux_plot, euv_flux_plot_94, euv_flux_plot_131, euv_flux_plot_193, euv_flux_plot_211, euv_flux_plot_304, euv_flux_plot_335,

    # euv_flux_plot.set_data(time_list, cropped_intensity_list)

    def update4(i):
        pass

    big_plot = FuncAnimation(fig2, update, frames=number_of_images, interval=21, blit=True)
    small_plot = FuncAnimation(fig3, update2, frames=number_of_images, interval=21, blit=True)
    graph_plot = FuncAnimation(fig4, update3, frames=number_of_images, interval=21, blit=True)
    # flare_point_player = FuncAnimation(fig3,point_update, frames=range(number_of_images), interval=10, blit=True)

    big_plot.running = True

    def pause():
        big_plot.event_source.stop()
        small_plot.event_source.stop()
        graph_plot.event_source.stop()
        big_plot.running = False

    def play():
        if big_plot.running == False:
            big_plot.event_source.start()
            small_plot.event_source.start()
            graph_plot.event_source.start()
            big_plot.running = True
        else:
            big_plot.event_source.stop()
            small_plot.event_source.stop()
            graph_plot.event_source.stop()
            big_plot.running = False

        # flare_point_player.event_source.start()

    def step():
        print("step used")
        global one_frame_state
        global current_frame_number
        print(current_frame_number)
        global graph_type_changed_bool
        global graph_type_regular_bool

        one_frame_state = True
        if big_plot.running == False:
            print("movie stopped")
            update(current_frame_number)
            update2(current_frame_number)
            # update3(current_frame_number)

            print("ax13 lines")

            print(ax13.lines)
            if len(ax13.lines) > 0:
                ax13.lines.clear()

            print(ax13.lines)
            print("ax12 lines")
            print(ax12.lines)
            if len(ax12.lines) > 0:
                ax12.lines.clear()
                # ax12.lines[0].remove()

            flare_x_center = original_solar_movies_dict[current_frame_number][0]
            flare_y_center = original_solar_movies_dict[current_frame_number][1]
            ax12.plot(y_center / 8, x_center / 8, 's', color='green', fillstyle='none', markersize=30)
            """
            if len(ax17.lines) > 0:
                ax17.lines.clear()
            """

            if key:
                pass
            else:
                ax12.plot(changing_solar_movies_dict["ix"], changing_solar_movies_dict["iy"], 's', color='red',
                          fillstyle='none', markersize=30)

            if (changing_solar_movies_dict["small_x"] < flare_x_center < changing_solar_movies_dict["big_x"]) and (
                    changing_solar_movies_dict["small_y"] < flare_y_center < changing_solar_movies_dict["big_y"]):
                ax13.plot(original_solar_movies_dict[current_frame_number][0],
                          original_solar_movies_dict[current_frame_number][1], 'wx', color='red', fillstyle='none',
                          markersize=15)

            canvas12.draw()
            canvas13.draw()
            # canvas14.draw()
            print("ax14 lines")
            print(ax14.lines)
            if len(ax14.lines) > 3:
                ax14.lines[3].remove()

            time = chosen_time_dict[current_frame_number]
            real_time = to_datetime_object(time, date_format)
            ax14.axvline(real_time, ls='-', color='black', lw=1, zorder=10)

            if graph_type_regular_bool:
                ax17.plot_date(time_list, cropped_intensity_list, marker='', linestyle='-', color='darkgreen',
                               label="Intensity")
            else:
                ax17.plot_date(time_list, normalized_cropped_intensity_list, marker='', linestyle='-',
                               color='darkgreen', label="Intensity")

            print("booleans")
            print(graph_type_regular_bool)
            print(graph_type_changed_bool)

            ax17.relim()
            ax17.autoscale_view()

            # fig4.tight_layout()
            if user_cropped_map_intensity_list:
                if graph_type_regular_bool:

                    ax17.plot_date(time_list, user_cropped_map_intensity_list, marker='', linestyle='-', color='red',
                                   label="User Intensity")
                else:
                    ax17.plot_date(time_list, changing_solar_movies_dict[
                        "user_selected_normalized_graph"], marker='', linestyle='-', color='red',
                                   label="User Intensity")
                # ax14.relim()
                # ax14.autoscale_view()

            if graph_type_changed_bool:
                print("redrew canvas")

                ax17.relim()
                ax17.autoscale_view()
                graph_type_changed_bool = False

            canvas14.draw()

            one_frame_state = False
            slider_val.set(current_frame_number)

    def one_forward():

        global current_frame_number
        print(current_frame_number)
        current_frame_number += 1
        if current_frame_number + 1 > number_of_images:
            current_frame_number = -1
        step()
        print(current_frame_number)

    def one_backward():
        global current_frame_number
        print("current frame number " + str(current_frame_number))
        current_frame_number -= 1
        print("current frame number " + str(current_frame_number))
        if current_frame_number < 0:
            current_frame_number = 0
        step()

    def movie_speed_slider_func(i):
        global animation_interval
        global animation_interval_changed
        slider_number = int(movie_speed_slider.get())
        animation_interval = int((-75 * slider_number) + 1075)
        animation_interval_changed = True
        print("speed slider")
        print(animation_interval)

        # past_input_dict[str(event_peak_time1)[:-7] + str(wavelength_value)] = [plot_1, submapsequence, no_spikes_list, hpc_max_x,
        #                                              hpc_max_y, exposure_time_list]

    global center_flare_label
    center_flare_label = Label(root, bg="white", text="Center of\nFlare Position", justify=CENTER)
    user_center_flare_label = Label(root, bg="white", text="Box Position", justify=CENTER)

    center_flare_label.place(x=810 + space, y=545)
    user_center_flare_label.place(x=810 + space, y=605)

    quit_button.place(x=825 + space, y=655)

    fig4.autofmt_xdate()

    rhessi_center_flare_label = Label(root, bg="white", text="RHESSI Center \nPosition", justify=CENTER)
    rhessi_user_center_flare_label = Label(root, bg="white", text="User RHESSI\nCenter Position", justify=CENTER)

    rhessi_labels_x = 810 + space
    rhessi_labels_y = 428
    rhessi_center_flare_label.place(x=rhessi_labels_x, y=rhessi_labels_y)
    rhessi_user_center_flare_label.place(x=rhessi_labels_x, y=rhessi_labels_y + 59)

    rhessi_flare_coordinates = Label(root, bg="white", textvariable=rhessi_flare_text)
    user_rhessi_flare_coordinates = Label(root, bg="white", textvariable=user_rhessi_flare_text)
    rhessi_flare_coordinates.place(x=rhessi_labels_x, y=rhessi_labels_y + 35)
    user_rhessi_flare_coordinates.place(x=rhessi_labels_x, y=rhessi_labels_y + 94)

    def flare_plot_norm_selection_for_button():
        global im, im2
        normalization = normalization_string.get()
        im = flare_plot_norm_selection(normalization, ax12, full_image_dict, colormap, [0, 512, 512, 0])
        im2 = flare_plot_norm_selection(normalization, ax13, cropped_image_list, colormap,
                                        [small_x, big_x, big_y, small_y])

    asinh_normalization_button.config(command=flare_plot_norm_selection_for_button)
    normal_normalization_button.config(command=flare_plot_norm_selection_for_button)

    def onclick(event):
        global user_clicked_map_bool
        # ani.stop()
        # ani3.stop()
        ax13.clear()

        user_cropped_map_intensity_list.clear()
        user_clicked_map_bool = True

        if (event.inaxes is not None):
            ix, iy = int(event.xdata) * 8, int(event.ydata) * 8
            if not iy == 100000000:
                key.clear()
                changing_solar_movies_dict["center"] = [ix, iy]
                cropped_image_list.clear()
                frame_counter = 0
                for frame in finalized_frame_list:
                    cropped_image_data = all_images_dict[frame][iy - 300:iy + 300, ix - 300: ix + 300]
                    cropped_image_list.append(cropped_image_data)
                    intensity_value = cropped_image_data.sum()
                    user_cropped_map_intensity_list.append(intensity_value)
                    frame_counter += 1

                # if len(ax17.lines) > 2:
                #   ax17.lines[-1].remove()

                small_x = ix - square_side_length
                small_y = iy - square_side_length
                big_x = ix + square_side_length
                big_y = iy + square_side_length

                helioprojective_small_x = x_pixel_to_helioprojective(small_x)
                helioprojective_small_y = y_pixel_to_helioprojective(small_y)
                rounded_small_x = roundup(helioprojective_small_x)
                rounded_small_y = rounddown(helioprojective_small_y)
                first_x = x_helioprojective_to_pixel(rounded_small_x, 0)
                second_x = x_helioprojective_to_pixel(rounded_small_x, 100)
                third_x = x_helioprojective_to_pixel(rounded_small_x, 200)
                first_y = y_helioprojective_to_pixel(rounded_small_y, 0)
                second_y = y_helioprojective_to_pixel(rounded_small_y, 100)
                third_y = y_helioprojective_to_pixel(rounded_small_y, 200)
                rect = patches.Rectangle((ix / 8 - square_side_length / 8, iy / 8 - square_side_length / 8),
                                         square_side_length / 4,
                                         square_side_length / 4, linewidth=1.5, edgecolor='g', facecolor='none')

                changing_solar_movies_dict["helioprojective_small_x"] = helioprojective_small_x
                changing_solar_movies_dict["helioprojective_small_y"] = helioprojective_small_y
                changing_solar_movies_dict["small_x"] = small_x
                changing_solar_movies_dict["small_y"] = small_y
                changing_solar_movies_dict["big_x"] = big_x
                changing_solar_movies_dict["big_y"] = big_y
                changing_solar_movies_dict["rounded_small_x"] = rounded_small_x
                changing_solar_movies_dict["rounded_small_y"] = rounded_small_y
                changing_solar_movies_dict["first_x"] = first_x
                changing_solar_movies_dict["second_x"] = second_x
                changing_solar_movies_dict["third_x"] = third_x
                changing_solar_movies_dict["first_y"] = first_y
                changing_solar_movies_dict["second_y"] = second_y
                changing_solar_movies_dict["third_y"] = third_y
                # changing_solar_movies_dict["rectangle_boarder"] = rect
                changing_solar_movies_dict["user_selected_center_pixel"] = ix, iy
                changing_solar_movies_dict["ix"] = ix / 8
                changing_solar_movies_dict["iy"] = iy / 8
                user_projective_x = x_pixel_to_helioprojective(ix)
                user_projective_y = y_pixel_to_helioprojective(iy)
                changing_solar_movies_dict[
                    "user_selected_center_helioprojective"] = user_projective_x, user_projective_y

                changing_solar_movies_dict["user_selected_center_heliographical"] = helioprojective_to_heliographic(
                    user_projective_x, user_projective_y, random_time)

                normalized_user_cropped_map_intensity_list = normalized_plot(user_cropped_map_intensity_list)
                changing_solar_movies_dict[
                    "user_selected_normalized_graph"] = normalized_user_cropped_map_intensity_list

                change_ax_13_axis()
                canvas13.draw_idle()

                normalization = normalization_string.get()
                global im2

                im2 = flare_plot_norm_selection(normalization, ax13, cropped_image_list, colormap,
                                                [small_x, big_x, big_y, small_y])

                h2, l2 = ax17.get_legend_handles_labels()
                ax14.legend(h1 + h2, l1 + l2, loc=2)
                ax17.relim()
                ax17.autoscale_view()
                big_plot.event_source.start()
                small_plot.event_source.start()
                graph_plot.event_source.start()
                big_plot.running = True
                cid1 = fig2.canvas.mpl_connect('button_press_event', onclick)




            else:
                print('Clicked ouside axes bounds but inside plot window')
        else:
            print('Clicked ouside axes bounds but inside plot window')

    def rhessi_onclick(event):
        print("onclick")
        ax14.clear()
        if event.inaxes is not None:
            ix, iy = int(event.xdata), int(event.ydata)
            print(ix, iy)
            if not iy == 100000000:

                rhessi_helioprojective_x = rhessi_pixel_to_helioprojective_x(rhessi_header_x, ix, rhessi_multiplier)
                rhessi_helioprojective_y = rhessi_pixel_to_helioprojective_y(rhessi_header_y, iy, rhessi_multiplier)
                rhessi_heliographical_coords = helioprojective_to_heliographic(rhessi_helioprojective_x,
                                                                               rhessi_helioprojective_y, random_time)

                rhessi_axis_dict["rhessi_helioprojective_x"] = rhessi_helioprojective_x
                rhessi_axis_dict["rhessi_helioprojective_y"] = rhessi_helioprojective_y
                rhessi_axis_dict["rhessi_heliographical_coords"] = rhessi_heliographical_coords

                coord_change_func()
                # ax16.lines[-1].remove()
                ax16.plot(ix, iy, 'wx', color='green', fillstyle='none', markersize=15)

                canvas13.draw()

                cid1 = fig1.canvas.mpl_connect('button_press_event', rhessi_onclick)

            else:
                print('Clicked ouside axes bounds but inside plot window')
        else:
            print('Clicked ouside axes bounds but inside plot window')

    def frame_slider_func(i):
        global current_frame_number
        global slider_used_func
        print("frame func used")
        slider_frame_number = int(slider_val.get())
        current_frame_number = slider_frame_number
        step()

    cid = fig2.canvas.mpl_connect('button_press_event', onclick)
    cid1 = fig1.canvas.mpl_connect('button_press_event', rhessi_onclick)

    play_button = Button(frame2, highlightbackground="white", bg="white", text="Play/Pause", command=play)

    one_forward_button = Button(frame2, highlightbackground="white", bg="white", text="One Foward", command=one_forward)

    one_backward_button = Button(frame2, highlightbackground="white", bg="white", text="One Backward",
                                 command=one_backward)

    slider_val = DoubleVar()

    movie_speed_slider = Scale(frame2, from_=1, to=15, orient=HORIZONTAL, highlightbackground="white", bg="white",
                               command=movie_speed_slider_func)
    frame_slider = Scale(frame2, from_=0, to=number_of_images, orient=HORIZONTAL, variable=slider_val,
                         highlightbackground="white", bg="white",
                         command=frame_slider_func, length=200)

    movie_speed_slider.set(14)
    movie_speed_slider.place(x=130, y=738)
    global current_frame_number
    slider_val.set(current_frame_number)
    frame_slider.place(x=420, y=385)

    play_button.place(x=260, y=750)
    one_forward_button.place(x=370, y=750)
    one_backward_button.place(x=480, y=750)
    fig1.tight_layout()

    fig4.tight_layout()

    GOES_graph_type_string = StringVar()
    GOES_graph_type_string.set("logarithmic")

    root.attributes("-topmost", True)


def get_time_range_function():
    global flare_selected
    flare_selected = False
    clicked_flare_list = []
    master_dict.clear()
    try:
        center_flare_label.place_forget()
    except NameError:
        pass

    try:
        flare_center_coordinates.place_forget()
    except NameError:
        pass
    try:
        frame2.destroy()
    except NameError:
        pass
    # frame1.destroy()
    create_frame2()
    # create_frame()

    for widgets in root.winfo_children():
        print("wdigeg")
        print(widgets)
        widgets.place_forget()

    global canvas11
    try:
        canvas11.get_tk_widget().pack_forget()
    except Exception:
        print("Error ignored, no figures yet")
    get_time_range_button_old.destroy()
    root.geometry("1200x1000")
    start_date = start_date_entry_box.get()
    end_date = end_date_entry_box.get()

    global tr
    if end_date == "":
        end_date = start_date
        tr = TimeRange([start_date, end_date])
        tr.extend(TimeDelta(0, format='sec'), TimeDelta(86400, format='sec'))
        end_date_string.set(str(tr.end)[:-13])

    else:

        tr = TimeRange([start_date, end_date])

    fig1, ax11 = plt.subplots(figsize=(8.3, 8))

    canvas11 = FigureCanvasTkAgg(fig1, frame1)
    # canvas11.get_tk_widget().config(highlightthickness = 0, borderwidth = 0, bd = 5)

    navigation_bar = tkagg.NavigationToolbar2Tk(canvas11, frame1)
    navigation_bar.place(x=1000, y=2000)

    try:
        canvas11.get_tk_widget().pack()
    except Exception:
        canvas11.get_tk_widget().place(x=0, y=0)

    client = hek.HEKClient()

    """last_flare_end_time = parse_time(flares_hek[len(flares_hek) - 1].get('event_endtime'))

    if last_flare_end_time > tr.end:
        tr = TimeRange([start_date, last_flare_end_time])
    """
    print("hiiiii")
    global goes

    fake_goes = pd.read_pickle("2011ALLDATA.pickle")['20110213':'20110216']
    goes = ts.TimeSeries(fake_goes, source='xrs', concatenate=True)

    global minpeaks
    global maxpeaks
    # goes = nontruncgoes.truncate(tr)
    print("goes")
    print(goes)

    goes_plot = goes.plot(axes=ax11)
    series = goes.to_dataframe()['xrsb']
    minpeaks, maxpeaks = findpeaks(series, DELTA=1e-7)

    minpeaks_copy = minpeaks
    print("maxpeaks")
    print(minpeaks_copy)
    print(maxpeaks)

    ax11.set_yscale('log')

    ax11.set_xlabel(' ')
    ax11.set_ylabel('$Watts/m^{-2}$')

    ax11.set_ylim(1.0E-9, 1.0E-3)

    X_label = Label(frame1, text="X", font=("DejaVu Sans", 50), bg="white")
    M_label = Label(frame1, text="M", font=("DejaVu Sans", 50), bg="white")
    C_label = Label(frame1, text="C", font=("DejaVu Sans", 50), bg="white")
    B_label = Label(frame1, text="B", font=("DejaVu Sans", 50), bg="white")
    A1_label = Label(frame1, text="A", font=("DejaVu Sans", 50), bg="white")
    A2_label = Label(frame1, text="A", font=("DejaVu Sans", 50), bg="white")

    def show_flare_class_labels():
        X_label.place(x=35, y=50)
        M_label.place(x=29, y=165)
        C_label.place(x=32, y=280)
        B_label.place(x=35, y=395)
        A1_label.place(x=35, y=510)
        A2_label.place(x=35, y=625)

    show_flare_class_labels()

    def hide_flare_class_labels():
        X_label.place_forget()
        M_label.place_forget()
        C_label.place_forget()
        B_label.place_forget()
        A1_label.place_forget()
        A2_label.place_forget()

    fig1.tight_layout()
    """
    if not flares_hek:
        no_flares_label = Label(root, text="No Solar Flares Found \n Please input a different time range")
        no_flares_label.place(x=670, y=200)
        return "No Solar Flares Found"
    """

    global flare_event_list
    global hover_flare_time_list
    flare_event_list = []
    hover_flare_time_list = []

    try:
        minimum_value = str(goes_minimum_value_entry_box.get())
    except ValueError:
        minimum_value = ""
    try:
        maximum_value = str(goes_maximum_value_entry_box.get())
    except ValueError:
        maximum_value = ""

    if minimum_value == "None":
        minimum_value = ""
    if maximum_value == "None":
        maximum_value = ""

    def flare_class_labeling(class1, peak_time, max_flux_value):
        ax11.annotate(class1, xy=(peak_time, maximum_point),
                      xytext=(peak_time, max_flux_value), ha='center')

    global flare_dropdown_list
    global maximum_point_list
    global peak_time_list
    global current_flare_start_and_end_times
    global flare_class_list
    global start_time_list
    flare_dropdown_list = []
    maximum_point_list = []
    peak_time_list = []
    start_time_list = []
    end_time_list = []
    current_flare_start_and_end_times = []
    flare_class_list = []

    def magnitude_measurement(flare):
        print("flare")
        print(flare)
        flare_value = flare[1]
        if flare_value > 1.1e-4:
            class_val = "X"
            magnitide = str(flare_value)[5] + "." + str(flare_value)[6]

            magnitide = magnitide.replace("0", "")

        elif flare_value > 1.0e-5:
            class_val = "M"
            magnitide = str(flare_value)[:3]
        elif flare_value > 1.0e-6:
            class_val = "C"
            magnitide = str(flare_value)[:3]
        elif flare_value > 1.0e-7:
            class_val = "B"
            magnitide = str(flare_value)[:3]
        else:
            return False
        flare_class_val = class_val + magnitide
        print(flare_class_val)

        flare_class_list.append(flare_class_val)
        return flare_class_val

    for event in range(0, len(maxpeaks)):
        event_peak_time = maxpeaks[event][0]

        try:
            if maxpeaks[event + 1][0] - event_peak_time < timedelta(minutes=5) or abs(
                    maxpeaks[event + 1][1] - maxpeaks[event][1]) < 5e-08:
                continue
            else:
                print(event_peak_time)
                print(maxpeaks[event + 1][0] - event_peak_time)

        except IndexError:
            continue

        flare_class = magnitude_measurement(maxpeaks[event])
        if not flare_class:

            pass
        else:
            if flare_class == "C6.9":
                print("C6.9 denied")
                print(maxpeaks[event + 1][0] - event_peak_time)
                print(abs(maxpeaks[event + 1][1] - maxpeaks[event][1]))
                print(maxpeaks[event + 1][1], maxpeaks[event][1])
                print("done")

            event_start_time = minpeaks[event][0]
            peak_time_list.append(event_peak_time)
            start_time_list.append(event_start_time)
            try:
                end_time_list.append(minpeaks[event + 1][0])
            except IndexError:
                pass
            print("event start time")
            print(minpeaks[event])

            flare_dropdown_list.append(flare_class + " " + str(event_peak_time)[:-10])
            flare_event_list.append(flare_class)
    print("flare dropdown list")
    print(flare_dropdown_list)

    flare_labels_to_be_deleted = []

    def between(peak_time_list_func, minpeaks_list):
        print("minpeaks_list2")
        print(minpeaks_list)
        flare_split_list = []

        n = 0
        for g in range(0, len(minpeaks_list)):
            if peak_time_list_func[0] > minpeaks_list[g][0]:
                n += 1
            else:
                for _ in range(0, n):
                    minpeaks_list.pop(0)
                break

        num = 0
        for i in range(1, len(peak_time_list_func)):
            subflare_list = []
            for j in range(0, len(minpeaks_list)):

                if peak_time_list_func[i - 1] < minpeaks_list[j][0] < peak_time_list_func[i]:
                    subflare_list.append(minpeaks_list[j])
                    num += 1
                else:

                    if subflare_list:
                        val = 100000000
                        for _ in range(0, num):
                            minpeaks_list.pop(0)
                        num = 0
                        minimum_flare = ""

                        for k in range(0, len(subflare_list)):

                            if subflare_list[k][1] < val:
                                val = subflare_list[k][1]
                                minimum_flare = subflare_list[k]

                        minimum_flare = list(str(minimum_flare[0]))
                        minimum_flare[10] = "T"
                        minimum_flare = "".join(minimum_flare)
                        minimum_flare = Time(minimum_flare, format='isot', scale='utc')
                        flare_split_list.append(minimum_flare)
                        break
                    else:
                        # delete flare dropdown menu value i-1

                        flare_dropdown_list.pop(i - 1)
                        flare_labels_to_be_deleted.append(i - 1)
                        break
        return flare_split_list

    global flare_splits
    flare_splits = between(peak_time_list, minpeaks_copy)
    print("minpeaks copy")
    print(flare_splits)
    print(minpeaks_copy)
    print(maxpeaks)

    def colored_time_range(color_name):
        try:
            ax11.axvspan(flare_splits[i].plot_date,
                         flare_splits[i + 1].plot_date,
                         alpha=0.2, color=color_name, label=flare_class, lw=0)
        except IndexError:
            ax11.axvspan(flare_splits[i].plot_date,
                         tr.end.plot_date,
                         alpha=0.2, color=color_name, label=flare_class, lw=0)

    # ax11.scatter(*zip(*flare_splits), color='red', label='min')

    flare_splits.append(tr.end)
    flare_splits.insert(0, tr.start)

    for i in range(0, len(flare_splits)):
        try:
            time_range = TimeRange(flare_splits[i], flare_splits[i + 1])
            hover_flare_time_list.append(time_range)
        except IndexError:
            pass

    initial_index = 0
    global axvspan_list
    axvspan_list = []
    global odd
    for i in range(0, len(flare_class_list)):
        flare_class = flare_class_list[i]
        z = max(flare_event_list)
        initial_index = flare_event_list.index(z)
        if odd == True:
            colored_time_range("gray")
            odd = False
        else:
            odd = True
            pass

        if flare_class == z:
            try:
                x = ax11.axvspan(flare_splits[i].plot_date,
                                 flare_splits[i + 1].plot_date,
                                 alpha=0.8, color='red', label=flare_class, lw=0)
            except IndexError:
                x = ax11.axvspan(flare_splits[i].plot_date,
                                 tr.end.plot_date,
                                 alpha=0.8, color='red', label=flare_class, lw=0)

            clicked_flare_list.append(x)

            current_flare_start_and_end_times.append(flare_splits[i])
            try:
                current_flare_start_and_end_times.append(flare_splits[i + 1])
            except IndexError:
                current_flare_start_and_end_times.append(tr.end)

    label_number = 0

    # different from the other one
    flares_dropdown_menu.clear()
    real_flare_class_list.clear()
    real_start_time_list = []
    real_end_time_list = []

    def flare_labelling(q):
        if minimum_value == "" and maximum_value == "":
            if flare_class > "C0":
                flare_class_labeling(flare_class, time_in_string, maximum_point)
                flares_for_rhessi_comparison.append(time_in_string)

        elif not minimum_value == "" and not maximum_value == "":
            if minimum_value < flare_class < maximum_value:
                flare_class_labeling(flare_class, time_in_string, maximum_point)
                flares_for_rhessi_comparison.append(time_in_string)
        elif minimum_value == "" and not maximum_value == "":
            if flare_class < maximum_value:
                flare_class_labeling(flare_class, time_in_string, maximum_point)
                flares_for_rhessi_comparison.append(time_in_string)
        elif not minimum_value == "" and maximum_value == "":
            if flare_class > minimum_value:
                flare_class_labeling(flare_class, time_in_string, maximum_point)
                flares_for_rhessi_comparison.append(time_in_string)
                flares_dropdown_menu.append(flare_dropdown_list[q])
                real_peak_time_list.append(peak_time_list[q])
                real_flare_class_list.append(flare_class_list[q])
                real_start_time_list.append(start_time_list[q])
                real_end_time_list.append(end_time_list[q])

    for i in range(0, len(flare_class_list)):
        event_peak_time3 = peak_time_list[i]
        event_beginning_time3 = start_time_list[i]

        flare_class = flare_class_list[i]
        # ax11.legend(loc=2)

        beginning_time_data = str(event_beginning_time3)[:-7]
        time_in_string = str(event_peak_time3)[:19]

        maximum_point = max(goes.data['xrsb'][time_in_string])
        maximum_point_list.append(maximum_point)

        flare_labelling(i)

    flare_dropdown_list_max_flare = flare_dropdown_list[initial_index]
    real_initial_index = flares_dropdown_menu.index(flare_dropdown_list_max_flare)
    global initial_menu_value
    initial_menu_value = StringVar()
    initial_menu_value.set(flares_dropdown_menu[real_initial_index])

    number = [initial_index]
    hover_status = ["yes"]
    print("flares for rhessi")

    print(flares_for_rhessi_comparison)

    print("hover flare time list")
    print(peak_time_list)

    def hover(event):
        if hover_status[0] == "yes":
            if (event.inaxes is not None):
                clicked_flare_time = mdates.num2date(event.xdata)
                for i in range(0, len(hover_flare_time_list)):
                    if clicked_flare_time in hover_flare_time_list[i]:
                        try:
                            if not number[0] == i:
                                try:
                                    axvspan_list[0].remove()
                                    axvspan_list.clear()
                                except IndexError:
                                    pass

                                ab = ax11.axvspan(flare_splits[i].plot_date,
                                                  flare_splits[i + 1].plot_date,
                                                  alpha=0.3, color='red', label=flare_class, lw=0)
                                axvspan_list.append(ab)
                                current_flare_start_and_end_times.clear()
                                current_flare_start_and_end_times.append(flare_splits[i])
                                current_flare_start_and_end_times.append(flare_splits[i + 1])
                                fig1.canvas.draw()
                                fig1.tight_layout()
                                number[0] = i
                            else:
                                pass
                        except TypeError:
                            print("first flare")
                        break
                    else:
                        pass
            else:
                print('Clicked ouside axes bounds but inside plot window')

    did = fig1.canvas.mpl_connect('motion_notify_event', hover)

    def click_flare_time(event):
        if (event.inaxes is not None):
            clicked_flare_time = mdates.num2date(event.xdata)
            for i in range(0, len(hover_flare_time_list)):
                if clicked_flare_time in hover_flare_time_list[i]:

                    nu = flares_dropdown_menu.index(flare_dropdown_list[i])
                    initial_menu_value.set(flares_dropdown_menu[nu])
                    auto_time_range3()
                    clicked_flare_list[0].remove()
                    clicked_flare_list.clear()

                    clicked_flare = ax11.axvspan(flare_splits[i].plot_date,
                                                 flare_splits[i + 1].plot_date,
                                                 alpha=0.8, color='red', lw=0)  # , label=flare_class)
                    clicked_flare_list.append(clicked_flare)
                    fig1.canvas.draw()
                    break
                else:
                    # print("no flares in time range")
                    pass
        else:
            print('Clicked ouside axes bounds but inside plot window')

    cid = fig1.canvas.mpl_connect('button_press_event', click_flare_time)

    flare_entry_label = Label(root, highlightbackground="white", bg="white", text="Flare:")
    wavelength_entry_label = Label(root, highlightbackground="white", bg="white", text="Wavelength:")

    flare_entry_label.place(x=690 + space, y=210)
    wavelength_entry_label.place(x=690 + space, y=240)

    global initial_wavelength
    initial_wavelength = StringVar(root)
    initial_wavelength.set(wavelength_list[2])
    wavelength_dropdown_menu = OptionMenu(root, initial_wavelength, *wavelength_list)
    wavelength_dropdown_menu.config(bg="white")
    wavelength_dropdown_menu.place(x=770 + space, y=240)
    global flares_event_dropdown_button
    flares_event_dropdown_button = Button(root, highlightbackground="white", bg="white", text="Confirm\nSelections",
                                          command=flares_event_dropdown_function)
    flares_event_dropdown_button.place(x=810 + space, y=400)
    start_date_entry_label.place(x=690 + space, y=10)
    end_date_entry_label.place(x=690 + space, y=55)

    start_date_entry_box.place(x=690 + space, y=30)
    end_date_entry_box.place(x=690 + space, y=75)
    get_time_range_button = Button(root, highlightbackground="white", bg="white", text="Get Flare Info",
                                   command=get_time_range_function)
    get_time_range_button.place(x=733 + space, y=168)
    goes_minimum_value_entry_label.config(text="Minimum Flare Class:")
    goes_maximum_value_entry_label.config(text="Maximum Flare Class:")
    goes_minimum_value_entry_label.place(x=690 + space, y=110)
    goes_maximum_value_entry_label.place(x=690 + space, y=140)
    goes_minimum_value_entry_box.config(width=4)
    goes_maximum_value_entry_box.config(width=4)
    goes_minimum_value_entry_box.place(x=830 + space, y=107)
    goes_maximum_value_entry_box.place(x=830 + space, y=138)
    quit_button.place(x=750 + space, y=460)
    panel.place(x=690 + space, y=685)
    panel1.place(x=800 + space, y=685)

    navigation_bar = tkagg.NavigationToolbar2Tk(canvas11, frame1)
    navigation_bar.place(x=1000, y=2000)

    def reset():
        fig1.canvas.toolbar.home()
        if GOES_graph_placeholder_list[0] == "linear":
            GOES_graph_placeholder_list.clear()
            GOES_graph_placeholder_list.append("linear")

            ax11.set_yscale('linear')

            max_value = max(maximum_point_list)
            # ax11.set_ylim(0, max_value+0.00001)
            ax11.yaxis.set_major_formatter(MathTextSciFormatter("%1.0e"))

            hide_flare_class_labels()
            ax11.tick_params(axis='y', labelsize=8)
            canvas11.draw_idle()

        else:
            show_flare_class_labels()

        # Fatal Python Error: PyEval_RestoreThread: NULL tstate

    def pan():
        fig1.canvas.toolbar.pan()
        hide_flare_class_labels()

    def zoom():
        fig1.canvas.toolbar.zoom()
        hide_flare_class_labels()

    reset_graph_button = Button(frame1, highlightbackground="white", bg="white", text="Reset", command=reset,
                                borderwidth=0, width=2)
    pan_graph_button = Button(frame1, highlightbackground="white", bg="white", text="Pan", command=pan, borderwidth=0,
                              width=2)
    zoom_graph_button = Button(frame1, highlightbackground="white", bg="white", text="Zoom", command=zoom,
                               borderwidth=0, width=2)

    reset_graph_button.place(x=30, y=763)
    pan_graph_button.place(x=80, y=763)
    zoom_graph_button.place(x=130, y=763)

    GOES_graph_type_string = StringVar()
    GOES_graph_type_string.set("logarithmic")

    def select_GOES_graph_type():
        GOES_graph_type = GOES_graph_type_string.get()
        if GOES_graph_type == "linear":
            GOES_graph_placeholder_list.clear()
            GOES_graph_placeholder_list.append("linear")
            ax11.set_yscale('linear')

            max_value = max(maximum_point_list)

            ax11.yaxis.set_major_formatter(MathTextSciFormatter("%1.0e"))

            hide_flare_class_labels()
            ax11.tick_params(axis='y', labelsize=8)
            ax11.set_ylim(-0.00001, max_value + 0.00001)
            canvas11.draw_idle()

        if GOES_graph_type == "logarithmic":
            GOES_graph_placeholder_list.clear()
            GOES_graph_placeholder_list.append("logarithmic")
            ax11.set_yscale('log')

            ax11.relim()
            ax11.set_ylim(1.0E-9, 1.0E-3)
            ax11.tick_params(axis='y', labelsize=10)
            canvas11.draw_idle()

            show_flare_class_labels()

    linear_button = Radiobutton(frame1, highlightbackground="white", bg="white", text='Linear',
                                variable=GOES_graph_type_string,
                                value='linear', command=select_GOES_graph_type, borderwidth=0)

    logarithmic_button = Radiobutton(frame1, highlightbackground="white", bg="white", text='Logarithmic',
                                     variable=GOES_graph_type_string,
                                     value='logarithmic', command=select_GOES_graph_type, borderwidth=0)

    linear_button.place(x=185, y=766)
    logarithmic_button.place(x=250, y=766)

    def add_to_manual_flare_list():
        flare_event_name = initial_menu_value.get()
        index = flare_dropdown_list.index(flare_event_name)
        time_range_start = flare_splits[index]
        time_range_end = flare_splits[index + 1]
        manual_selected_time_range.clear()
        manual_selected_time_range.append(time_range_start)
        manual_selected_time_range.append(time_range_end)

    start_end_line = []

    time_range_start_string = StringVar()
    time_range_end_string = StringVar()
    range_start_time = StringVar()
    range_end_time = StringVar()
    range_state_str = StringVar()

    time_range_start_label = Label(root, bg="white", textvariable=time_range_start_string)
    time_range_end_label = Label(root, bg="white", textvariable=time_range_end_string)
    range_start_time_label = Label(root, bg="white", textvariable=range_start_time)
    range_end_time_label = Label(root, bg="white", textvariable=range_end_time)
    time_range_start_label.place(x=768 + space, y=305)
    time_range_end_label.place(x=762 + space, y=325)
    range_start_time_label.place(x=690 + space, y=305)
    range_end_time_label.place(x=690 + space, y=325)

    def write_flare_range_time_values():
        time_range_start = manual_selected_time_range[0]
        time_range_end = manual_selected_time_range[1]
        time_range_start_string.set(str(time_range_start)[:16].replace("T", ' '))
        time_range_end_string.set(str(time_range_end)[:16].replace("T", ' '))
        range_start_time.set("Range Start:")
        range_end_time.set("Range End:")
        # manual_selected_time_range.clear()

    add_to_manual_flare_list()
    write_flare_range_time_values()

    def select_flare_range():
        hover_status[0] = "no"
        start_end_line.clear()
        if len(ax11.lines) > 2:
            ax11.lines[-1].remove()
            ax11.lines[-1].remove()
        start_time = current_flare_start_and_end_times[0]
        end_time = current_flare_start_and_end_times[1]

        start_line = draggableline(ax11, start_time.plot_date)
        end_line = draggableline(ax11, end_time.plot_date)
        start_end_line.append(start_line)
        start_end_line.append(end_line)
        auto_verification.clear()
        pass

    def finish_select_flare_range():
        hover_status[0] = "yes"
        start_line = start_end_line[0]
        end_line = start_end_line[1]
        manual_selected_time_range.clear()
        time_range_start = mdates.num2date(start_line.XorY)
        time_range_end = mdates.num2date(end_line.XorY)
        manual_selected_time_range.append(str(time_range_start)[:-9].replace(" ", 'T'))
        manual_selected_time_range.append(str(time_range_end)[:-9].replace(" ", 'T'))

        if len(ax11.lines) > 2:
            ax11.lines[-1].remove()
            ax11.lines[-1].remove()
        pass
        canvas11.draw_idle()
        write_flare_range_time_values()
        auto_verification.clear()
        range_state_str.set("Manual Selection")

    global number_of_frames_box
    number_of_frames_string = StringVar()

    number_of_frames_box = Entry(root, highlightbackground="white", bg="white", text="25", width=5)
    number_of_frames_text = Label(root, highlightbackground="white", bg="white", text="Number of Frames:")
    number_of_frames_box.delete(0, END)
    number_of_frames_box.insert(0, "25")

    def auto_time_range():
        hover_status[0] = "yes"
        if len(ax11.lines) > 2:
            ax11.lines[-1].remove()
            ax11.lines[-1].remove()
        canvas11.draw_idle()
        number_of_frames_box.delete(0, END)
        number_of_frames_box.insert(0, "25")
        auto_verification.append("Yes")
        add_to_manual_flare_list()
        write_flare_range_time_values()
        range_state_str.set("Auto Selection")

    def auto_time_range2(value):
        hover_status[0] = "yes"
        if len(ax11.lines) > 2:
            ax11.lines[-1].remove()
            ax11.lines[-1].remove()
        canvas11.draw_idle()
        number_of_frames_box.delete(0, END)
        number_of_frames_box.insert(0, "25")
        auto_verification.append("Yes")
        add_to_manual_flare_list()
        write_flare_range_time_values()
        range_state_str.set("Auto Selection")

    def auto_time_range3():
        if hover_status[0] == "yes":
            if len(ax11.lines) > 2:
                ax11.lines[-1].remove()
                ax11.lines[-1].remove()
            # canvas11.draw_idle()
            auto_verification.append("Yes")
            add_to_manual_flare_list()
            write_flare_range_time_values()
            range_state_str.set("Auto Selection")

    auto_time_range()

    flare_dropdown_menu = OptionMenu(root, initial_menu_value, *flares_dropdown_menu, command=auto_time_range2)
    flare_dropdown_menu.config(bg="white")
    flare_dropdown_menu.config(width=16)
    flare_dropdown_menu.place(x=730 + space, y=210)

    select_flare_range_variable = StringVar()
    select_flare_range_variable.set("select")

    select_flare_range_button = Button(frame1, highlightbackground="white", bg="white", text='Manual Select Range',
                                       command=select_flare_range, borderwidth=0, width=12)
    finish_select_flare_range_button = Button(frame1, highlightbackground="white", bg="white", text='Finish Selection',
                                              command=finish_select_flare_range, borderwidth=0, width=8)
    new_selection_button = Button(frame1, highlightbackground="white", bg="white", text="Auto Select",
                                  command=auto_time_range, borderwidth=0, width=6)

    select_flare_range_button.place(x=347, y=763)
    finish_select_flare_range_button.place(x=488, y=763)
    new_selection_button.place(x=595, y=763)
    number_of_frames_text.place(x=690 + space, y=270)
    number_of_frames_box.place(x=815 + space, y=267)

    range_state_label = Label(root, highlightbackground="white", bg="white", text="Range State:")
    range_state = Label(root, highlightbackground="white", bg="white", textvariable=range_state_str)
    range_state_str.set("Auto Selection")

    range_state_label.place(x=690 + space, y=345)
    range_state.place(x=775 + space, y=345)

    global change_frame_func

    def change_frame_func():
        global GOES_frame
        print("changed frame")
        if GOES_frame == False:
            change_frame_button['text'] = 'Solar Movies'
            frame1.tkraise()
            GOES_frame = True
        else:
            change_frame_button['text'] = 'GOES Plot'
            frame2.tkraise()
            GOES_frame = False

    frame1.tkraise(frame2)

    global change_frame_button
    change_frame_button = Button(root, highlightbackground="white", bg="white", text="Solar Movies",
                                 command=change_frame_func)

    # change_frame_button.place(x=752, y=415)
    global normalization_string
    normalization_string = StringVar()
    normalization_string.set("asinh")

    global normal_normalization_button
    global asinh_normalization_button

    normal_normalization_button = Radiobutton(root, highlightbackground="white", bg="white", text='Normal',
                                              variable=normalization_string,
                                              value="None", borderwidth=0)

    asinh_normalization_button = Radiobutton(root, highlightbackground="white", bg="white", text='AsinhStretch',
                                             variable=normalization_string,
                                             value="asinh", borderwidth=0)

    normal_normalization_button.config(command=do_nothing)
    asinh_normalization_button.config(command=do_nothing)

    image_scaling_label = Label(root, bg="white", text="Image Scaling")
    normalization_placement = 380
    normal_normalization_button.place(x=690 + space, y=normalization_placement + 20)
    asinh_normalization_button.place(x=690 + space, y=normalization_placement + 40)

    image_scaling_label.place(x=703 + space, y=normalization_placement - 5)

    # _______________________________________-

    smaller_solar_dict.clear()

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


get_time_range_button_old = Button(root, highlightbackground="white", bg="white", text="Get Solar Flare Information",
                                   command=get_time_range_function)
get_time_range_button_old.place(x=30, y=195)


def quit():
    sys.exit(0)


quit_button = Button(root, text="Quit", bg="white", highlightbackground="white", command=quit, borderwidth=0)
quit_button.place(x=265, y=195)

root.mainloop()