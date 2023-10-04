from __future__ import absolute_import, division, print_function

from UiMainWindowFinal import Ui_MainWindow
from double_range_slider import QRangeSlider
import matplotlib
from scipy import ndimage
import numpy as np
import cv2
import pfsspy
import matplotlib.patches as patches
from astropy.wcs import WCS
from PIL import Image
from numpy import asarray
import matplotlib.pyplot as plt
import bisect
import sunpy.data.sample
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
matplotlib.use('QT5Agg')
from sunpy.visualization.colormaps import color_tables as ct
from numpy.lib.stride_tricks import as_strided
import matplotlib.ticker as mticker
from sunpy.map import Map
from bisect import bisect_left

import math
from astropy.visualization import AsinhStretch

from astropy.visualization.mpl_normalize import ImageNormalize

import sunpy.timeseries as ts

from datetime import datetime, timedelta
from astropy.table import Table
import sys
import matplotlib

matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import pandas as pd
import time as lmao
import copy
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy.wcs
from reproject import reproject_interp
import sunpy.map
from pffspy_util import car_to_cea
from sunpy.coordinates import frames
import pyqtgraph as pg
from PyQt5.QtGui import *
from PyQt5.QtCore import *
import pyqtgraph.ptime as ptime
from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.QtWidgets import *
from astropy.time import Time
from aiapy.calibrate import degradation
import os

flare_selected = False
one_frame_state = False
graph_type_regular_bool = True
graph_type_changed_bool = False
slider_changed = False
animation_interval_changed = False
user_clicked_map_bool = False
coord_type_changed_bool = True
playing = True

os.environ['QT_MAC_WANTS_LAYER'] = '1'

translate = QtCore.QCoreApplication.translate

changing_solar_movies_dict = {}

def roundup(x):
    if x > 0:
        return int(math.ceil(x / 100.0)) * 100
    else:
        return (int(math.ceil(x / 100.0)) * 100) - 100


def rounddown(x):
    if x > 0:
        return (int(math.ceil(x / 100.0)) * 100) - 100
    else:
        return int(math.ceil(x / 100.0)) * 100


def x_helioprojective_to_pixel(value, additional_value):
    return ((value + additional_value) / pixel_arcsec_x_list_placeholder[0]) + sun_center_x_list_placeholder[0]


def y_helioprojective_to_pixel(value, additional_value):
    return ((value + additional_value) / pixel_arcsec_y_list_placeholder[0]) + sun_center_y_list_placeholder[0]


def x_pixel_to_helioprojective(value, sun_center,cdelt):
    return (value - sun_center) * cdelt


def y_pixel_to_helioprojective(value, sun_center,
                               cdelt):
    return (value - sun_center) * cdelt


def RHESSI_pixel_to_helioprojective(value, header_center, multiplier):
    return (value - header_center) * multiplier


def RHESSI_helioprojective_to_pixel(value, header_center, multiplier, additional_value):
    return ((value + additional_value) / multiplier) + header_center


def helioprojective_to_heliographic_simplified(x, y, sunpy_map):
    while True:
        try:
            print(x,y)
            c = SkyCoord(x * u.arcsec, y * u.arcsec, frame=sunpy_map.coordinate_frame)

            heliographic_coord = (c.transform_to(frames.HeliographicStonyhurst))
            _ = float(str(heliographic_coord.lat).split('d')[0])
            EW_lon, NS_lat = heliographic_coord.to_string('decimal').split(' ')
            EW_lon = float(EW_lon)
            NS_lat = float(NS_lat)
            break

        except ValueError:
            if abs(x) > abs(y):
                x = x - 1 if x > 0 else x + 1
                y = y - abs(y / x) if y > 0 else y + abs(y / x)
            else:
                y = y - 1 if y > 0 else y + 1
                x = x - abs(x / y) if x > 0 else x + abs(x / y)

    return NS_lat, EW_lon, heliographic_coord

def simplify_heliographic_coords(helio_coord):
    if helio_coord == "N/A":
        return helio_coord
    if 'W' in helio_coord:
        first,second = helio_coord.split('W')
        return first[:5] + "W"+ second[:4]
    else:
        first, second = helio_coord.split('E')
        return first[:5] + "E"+ second[:4]

def helioprojective_to_heliographic(x, y, sunpy_map):
    NS_lat, EW_lon, _ = helioprojective_to_heliographic_simplified(x, y, sunpy_map)
    if NS_lat > 0:
        first_part = "N" + str(NS_lat)
    else:
        first_part = "S" + str(abs(NS_lat))

    if EW_lon > 0:
        second_part = "W" + str(EW_lon)
    else:
        second_part = "E" + str(abs(EW_lon))
    return first_part + second_part


def create_difference_and_ratio_lists(image_data_list):
    run_diff_list = []
    base_diff_list = []
    run_ratio_list = []
    base_ratio_list = []
    zero_image = image_data_list[0]
    zero_image = np.clip(strided_rescale(zero_image, 8),1,999999999999)
    for i in range(1, len(image_data_list)):
        start = lmao.time()
        first_image = image_data_list[i]
        second_image = image_data_list[i - 1]
        first_image = np.clip(strided_rescale(first_image, 8),1,999999999999)
        second_image = np.clip(strided_rescale(second_image, 8),1,999999999999)

        run_diff_image = np.subtract(first_image, second_image)
        base_diff_image = np.subtract(first_image, zero_image)
        run_ratio_image = np.divide(first_image, second_image)
        base_ratio_image = np.divide(first_image, second_image)

        run_diff_list.append(run_diff_image)
        base_diff_list.append(base_diff_image)
        run_ratio_list.append(run_ratio_image)
        base_ratio_list.append(base_ratio_image)

    return run_diff_list, base_diff_list, run_ratio_list, base_ratio_list

def make_HMI_image_list(image_list, min_clip, max_clip):
    new_image_list = []
    for image in image_list:
        image = np.clip(image, min_clip, max_clip)
        new_image_list.append(image)
    return new_image_list

def resize(a, wanted_size):
    b = np.zeros((wanted_size, wanted_size))

    for i in range(wanted_size):
        for j in range(wanted_size):
            idx1 = int(i * len(a) / wanted_size)
            idx2 = int(j * len(a) / wanted_size)
            b[i][j] = a[idx1][idx2]

    return b

def HMI_weighted_centering_func(array1, string, max_clip_val):
    new_array = array1
    pd.to_pickle(new_array,'new_array.pickle')
    if string == "greater":
        new_array[new_array < max_clip_val] = 0

    if string == "less":
        new_array = new_array * -1

        new_array[new_array < max_clip_val] = 0
    print(new_array)
    print("new array")
    y_coord, x_coord = ndimage.measurements.center_of_mass(new_array)


    return x_coord, y_coord

def clickBox(self, state):
    if state == QtCore.Qt.Checked:
        print('Checked')
    else:
        print('Unchecked')


def convert_to_Qpoint_list(list1, divide_factor):
    new_list = []
    for path in list1:
        my_new_list = [QtCore.QPoint(coord[0] / divide_factor, coord[1] / divide_factor) for coord in path]
        new_list.append(my_new_list)
    return new_list


class MyLabel(QLabel):
    def __init__(self):
        super(MyLabel, self).__init__()

    def paintEvent(self, event):
        super(MyLabel, self).paintEvent(event)
        pos = QPoint(50, 50)
        painter = QPainter(self)
        painter.drawText(pos, 'hello,world')
        painter.setPen(QColor(255, 255, 255))


class RectItem(QtWidgets.QGraphicsRectItem):
    def paint(self, painter, option, widget=None):
        super(RectItem, self).paint(painter, option, widget)
        painter.save()
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.setBrush(QtCore.Qt.red)
        painter.setPen(QtCore.Qt.green)
        # painter.drawEllipse(option.rect)
        painter.restore()


class PolygonItem(QtWidgets.QGraphicsPolygonItem):
    def paint(self, painter, option, widget=None):
        super(PolygonItem, self).paint(painter, option, widget)
        painter.save()
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.setBrush(QtCore.Qt.red)
        painter.setPen(QtCore.Qt.green)
        # painter.drawEllipse(option.rect)
        painter.restore()


class LineItem(QtWidgets.QGraphicsLineItem):
    def paint(self, painter, option, widget=None):
        super(LineItem, self).paint(painter, option, widget)
        painter.save()
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.setBrush(QtCore.Qt.red)
        painter.setPen(QtCore.Qt.green)
        # painter.drawEllipse(option.rect)
        painter.restore()


class MyArrowItem(pg.ArrowItem):
    def paint(self, p, *args):
        p.translate(-self.boundingRect().center())
        pg.ArrowItem.paint(self, p, *args)


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


def create_flux_curve(time_table, total_intensity):
    tbl_meta = {'t_key': 't_value'}
    table = Table([time_table, total_intensity], names=['time', 'intensity'], meta=tbl_meta)
    table.add_index('time')

    ts_table = ts.TimeSeries(table)

    return ts_table

def find_brightest_pixel(image_list, pixel_dict, add_x, add_y, sun_center_x,sun_center_y,pixel_arcsec_x,pixel_arcsec_y, aia_map_object1):
    for i in range(0, len(image_list)):
        image_data = image_list[i]
        pixel_pos = np.argwhere(image_data == image_data.max()) * u.pixel
        point = pixel_pos[int(len(pixel_pos) / 2)]
        x = int(float(str(point[0])[:-4])) + add_x
        y = (int(float(str(point[1])[:-4])) + add_y)

        helioprojective_x = x_pixel_to_helioprojective(x, sun_center_x, pixel_arcsec_x)
        helioprojective_y = y_pixel_to_helioprojective(y, sun_center_y, pixel_arcsec_y)

        if helioprojective_x > 1100 or helioprojective_x < -1100 or helioprojective_y > 1100 or helioprojective_y < -1100:
            helioprojective_x = 0
            helioprojective_y = 0

        heliographic_coords = helioprojective_to_heliographic(helioprojective_x, helioprojective_y,
                                                              aia_map_object1)


        pixel_dict[i] = {}
        pixel_dict[i]['pixel'] = x, y
        pixel_dict[i]['cropped_pixel'] = x - add_x, y - add_y
        pixel_dict[i]['helioprojective'] = helioprojective_x, helioprojective_y
        pixel_dict[i]['heliographic'] = heliographic_coords


def convert_time_list_to_time_string_list(time_list):
    time_string_list = []
    for date1 in time_list:
        time_string_list.append(date1.strftime("%Y-%m-%d %H:%M:%S"))
    return time_string_list

def create_polygon(linepath, color):
    polygon = PolygonItem(QtGui.QPolygonF(linepath))
    polygon.setPen(color)
    polygon.setBrush(color)
    polygon.setOpacity(0.3)
    return polygon

def strided_rescale(g, bin_fac):
    strided = as_strided(g,
                         shape=(g.shape[0] // bin_fac, g.shape[1] // bin_fac, bin_fac, bin_fac),
                         strides=((g.strides[0] * bin_fac, g.strides[1] * bin_fac) + g.strides))
    return strided.mean(axis=-1).mean(axis=-1)  # order is NOT important! See notes..

def asinh_stretch_images(sdoaia, clip_val_first, clip_val_second, *argv):
    final_list = []
    for ratio_list in argv:

        new_list = []
        for image in ratio_list:
            image = np.clip(image, clip_val_first, clip_val_second)
            image = ImageNormalize(stretch=AsinhStretch(
                0.01)).__call__(image)
            image = (np.uint8(sdoaia(image) * 255))
            new_list.append(image)
        if len(argv) == 1:
            return new_list
        final_list.append(new_list)
    return final_list

def normalized_plot(plot_list):
    new_list = []

    max_value = max(plot_list)
    for intensity in plot_list:
        new_val = intensity / max_value
        new_list.append(new_val)
    return new_list

def to_datetime_object(date_string, date_format):
    s = datetime.strptime(str(date_string), date_format)
    return s

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


def HMI_to_CEA_func(m_hmi, y_resolution, x_resolution):
    """
    sun_center_x_list_placeholder = []
    sun_center_y_list_placeholder = []
    pixel_arcsec_x_list_placeholder = []
    pixel_arcsec_y_list_placeholder = []
    """

    beginning = lmao.time()
    shape_out = [y_resolution, x_resolution]
    frame_out = SkyCoord(0, 0, unit=u.deg, rsun=m_hmi.coordinate_frame.rsun, frame="heliographic_stonyhurst",
                         obstime=m_hmi.date)

    header = sunpy.map.make_fitswcs_header(
        shape_out,
        frame_out,
        scale=[180 / shape_out[0], 360 / shape_out[1]] * u.deg / u.pix,
        projection_code="CAR",
    )

    out_wcs = astropy.wcs.WCS(header)
    array, _ = reproject_interp(m_hmi, out_wcs, shape_out=shape_out)

    array = np.where(np.isnan(array), 0, array)
    m_hmi_cea = car_to_cea(sunpy.map.Map(array, header))

    m_hmi_cea.meta['TELESCOP'] = m_hmi.meta['TELESCOP']
    m_hmi_cea.meta['CONTENT'] = 'Carrington Synoptic Chart Of Br Field'
    m_hmi_cea.meta['T_OBS'] = m_hmi_cea.meta.pop(
        'DATE-OBS')  # This is because of a bug where the date accidentally returns None if it is in the date-obs key

    m_hmi_cea = sunpy.map.Map(
        m_hmi_cea.data,
        m_hmi_cea.meta,
    )

    return m_hmi_cea, out_wcs


def get_peak_HMI_data(m_hmi, x_center, y_center):
    def unlock_heliographic_coords(x, y):
        heliographic_coord = WCS(m_hmi_cea.meta).pixel_to_world(x, y)


        NS_lat = float(str(heliographic_coord.lat).replace('.', '').replace('d', '.').replace('m', '').replace('s', ''))
        EW_lon = float(str(heliographic_coord.lon).replace('.', '').replace('d', '.').replace('m', '').replace('s', ''))

        return NS_lat, EW_lon

    m_hmi_cea, out_wcs = HMI_to_CEA_func(m_hmi, 2048, 4096)
    m_hmi_cea_copy = copy.deepcopy(m_hmi_cea)

    full_data = (m_hmi_cea.data)
    full_data2 = copy.deepcopy(full_data)


    x_center_helioprojective = x_pixel_to_helioprojective(x_center, sun_center_x_list_placeholder[0], pixel_arcsec_x_list_placeholder[0])
    y_center_helioprojective = y_pixel_to_helioprojective(y_center,  sun_center_y_list_placeholder[0], pixel_arcsec_y_list_placeholder[0])

    y_center_heliographical_coord, x_center_heliographical_coord, heliographic_skycoord = helioprojective_to_heliographic_simplified(
        int(x_center_helioprojective), int(y_center_helioprojective), m_hmi)
    b = m_hmi_cea_copy.world_to_pixel(heliographic_skycoord)

    x_coord = int(float(str(b[0])[:-3]))
    y_coord = int(float(str(b[1])[:-3]))

    bruh = 300

    placeholder_array = np.zeros(shape=(500, 500))
    positive_array = full_data[y_coord - bruh:y_coord + bruh,
                     x_coord - bruh:x_coord + bruh]

    negative_array = copy.deepcopy(positive_array)

    x_positive, y_positive = HMI_weighted_centering_func(positive_array, 'greater', 50)
    x_negative, y_negative= HMI_weighted_centering_func(negative_array, 'less', 50)


    x_positive = x_positive + x_coord - 300
    x_negative = x_negative + x_coord - 300
    y_negative = y_negative + y_coord - 300
    y_positive = y_positive + y_coord - 300

    x_negative_heliographic, y_negative_heliographic = unlock_heliographic_coords(x_negative, y_negative)
    x_positive_heliographic, y_positive_heliographic = unlock_heliographic_coords(x_positive, y_positive)
    hmi_heliographic_x = (x_positive_heliographic + x_negative_heliographic) / 2
    hmi_heliographic_y = (y_positive_heliographic + y_negative_heliographic) / 2
    # x_positive_helioprojective, y_positive_helioprojective = heliographic_to_helioprojective(hmi_heliographic_x,
    #                                                                                        hmi_heliographic_y,
    #                                                                                       '2011-02-15')

    hmi_center_x_pixel = (x_negative + x_positive) / 2
    hmi_center_y_pixel = (y_negative + y_positive) / 2

    return m_hmi_cea, x_negative, x_positive, y_negative, y_positive, hmi_center_x_pixel, hmi_center_y_pixel, full_data2, m_hmi_cea_copy, out_wcs, x_center_heliographical_coord, y_center_heliographical_coord


wavelength_list =  [1600, 94, 171, 211, 304]

def create_image_list_dict(wavelength_dict, full_raw_image_list, wavelength, peak_x, peak_y, max_clip_val, date_list,
                           exposure_list):
    colormap = ct.aia_color_table(wavelength * u.angstrom)
    intensity_list = []
    raw_image_data_list = []
    raw_cropped_image_list = []

    bruh = 300
    t = Time(date_list[0], scale='utc')
    correction_factor = float(degradation(wavelength * u.angstrom, t)[0])
    for i in range(0, len(full_raw_image_list)):
        array = full_raw_image_list[i]
        cropped_image_data = array[peak_x - bruh:peak_x + bruh,
                             peak_y - bruh:peak_y + bruh]
        intensity_value = np.divide(cropped_image_data.sum(), exposure_list[i])
        intensity_list.append(intensity_value/correction_factor)
        raw_cropped_image_list.append(cropped_image_data)

        low_resolution_image_data = strided_rescale(array, 8)
        raw_image_data_list.append(low_resolution_image_data)

    image_data_list = asinh_stretch_images(colormap, 0, max_clip_val, raw_image_data_list)
    alt_cropped_image_list = asinh_stretch_images(colormap, 0, max_clip_val, raw_cropped_image_list)
    wavelength_dict[wavelength]["max_clip_val"] = max_clip_val

    raw_run_diff_list, raw_base_diff_list, raw_run_ratio_list, raw_base_ratio_list = create_difference_and_ratio_lists(
        full_raw_image_list)

    run_diff_list, base_diff_list = asinh_stretch_images(
        colormap, 0, 1500, raw_run_diff_list, raw_base_diff_list)

    run_ratio_list, base_ratio_list = asinh_stretch_images(
        colormap, 0, 10, raw_run_ratio_list, raw_base_ratio_list)

    time_string_list = convert_time_list_to_time_string_list(date_list)

    wavelength_dict[wavelength]['raw_run_diff_list'] = raw_run_diff_list
    wavelength_dict[wavelength]['raw_base_diff_list'] = raw_base_diff_list
    wavelength_dict[wavelength]['raw_run_ratio_list'] = raw_run_ratio_list
    wavelength_dict[wavelength]['raw_base_ratio_list'] = raw_base_ratio_list
    wavelength_dict[wavelength]['run_diff_list'] = run_diff_list
    wavelength_dict[wavelength]['base_diff_list'] = base_diff_list
    wavelength_dict[wavelength]['run_ratio_list'] = run_ratio_list
    wavelength_dict[wavelength]['base_ratio_list'] = base_ratio_list
    wavelength_dict[wavelength]["intensity_list"] = intensity_list
    wavelength_dict[wavelength]["normalized_intensity_list"] = normalized_plot(intensity_list)
    wavelength_dict[wavelength]["full_raw_image_data"] = full_raw_image_list
    wavelength_dict[wavelength]["raw_image_data"] = raw_image_data_list
    wavelength_dict[wavelength]["image_data"] = image_data_list
    wavelength_dict[wavelength]["date_list"] = date_list
    wavelength_dict[wavelength]["time_string_list"] = time_string_list
    wavelength_dict[wavelength]['raw_cropped_image_data'] = raw_cropped_image_list
    wavelength_dict[wavelength]["cropped_image_data"] = alt_cropped_image_list

#file_name = pd.read_pickle('pickle_file_name.pickle')
file_name = '/Users/jamisenma/Library/Application Support/JetBrains/PyCharmCE2021.1/scratches/Flare_Pickle_Files/M2.5 2011-06-07 06:41.pickle'
brightest_pixel_dict, peak_x, peak_y, RHESSI_center_x, RHESSI_center_y, all_wavelength_dict, non_wavelength_list, AR_dict, special_beginning_time_str, special_end_time_str, hmi_map, aia_map_object, flare_metadata = pd.read_pickle(
   file_name)  # add aia map variable later

sun_center_x_list_placeholder = [all_wavelength_dict[171]['sun_center_x']]
sun_center_y_list_placeholder = [all_wavelength_dict[171]['sun_center_y']]
pixel_arcsec_x_list_placeholder = [all_wavelength_dict[171]['pixel_arcsec_x']]
pixel_arcsec_y_list_placeholder = [all_wavelength_dict[171]['pixel_arcsec_y']]
sun_radius_list_placeholder = [all_wavelength_dict[171]['sun_radius']]

HMI_object = get_peak_HMI_data(hmi_map, peak_x, peak_y)



for wavelength in wavelength_list:
    full_raw_image_list = all_wavelength_dict[wavelength]['full_raw_image_data']
    datamin = all_wavelength_dict[wavelength]['datamin']
    date_list = all_wavelength_dict[wavelength]['date_list']
    exposure_list = all_wavelength_dict[wavelength]['exposure_list']
    create_image_list_dict(all_wavelength_dict, full_raw_image_list, wavelength, peak_x, peak_y, datamin, date_list,
                       exposure_list)

find_brightest_pixel(all_wavelength_dict[171]['raw_cropped_image_data'], brightest_pixel_dict, peak_x - 300,
                     peak_y - 300, sun_center_x_list_placeholder[0],sun_center_y_list_placeholder[0], pixel_arcsec_x_list_placeholder[0],pixel_arcsec_y_list_placeholder[0], aia_map_object)

event_beginning_time1 = flare_metadata['GOES Start Time']
event_peak_time1 = flare_metadata['GOES Peak Time']
event_end_time1 = flare_metadata['GOES End Time']

all_wavelength_dict['HMI']['max_clip_val'] = 50
all_wavelength_dict['RHESSI']['max_clip_val'] = 100



class ImageView(pg.ImageView):

    # constructor which inherit original
    # ImageView
    def __init__(self, *args, **kwargs):
        pg.ImageView.__init__(self, *args, **kwargs)


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=2, height=4, dpi=100):
        self.fig = Figure(figsize=(3, 2), dpi=dpi)
        super(MplCanvas, self).__init__(self.fig)


class AIA_Image_Widgets():
    def __init__(self, win, window, wavelength, size):
        self.win = win
        pg.setConfigOptions(antialias=True)
        self.plot = pg.PlotItem()
        self.window = window
        self.time_label = pg.TextItem('helo', **{'color': '#ffffff'})
        self.wavelength_label = pg.TextItem('helo', **{'color': '#ffffff'})
        self.wavelength = wavelength
        self.time_string_list = all_wavelength_dict[self.wavelength]['time_string_list']
        # self.max_clip_val = all_wavelength_dict[self.wavelength]['max_clip_val']
        self.sun_center_x = all_wavelength_dict[self.wavelength]['sun_center_x']
        self.sun_center_y = all_wavelength_dict[self.wavelength]['sun_center_y']
        self.pixel_arcsec_x = all_wavelength_dict[self.wavelength]['pixel_arcsec_x']
        self.pixel_arcsec_y = all_wavelength_dict[self.wavelength]['pixel_arcsec_y']
        self.size = size
        self.state = 'Normal'
        self.max_clip_val = int(all_wavelength_dict[self.wavelength]['max_clip_val'])
        if (wavelength == 'HMI' and size == 'cropped') or (wavelength == 1600 and size == 'cropped'):
            pass

        elif wavelength == 'HMI' and size == 'full':
            self.state = 'Normal'
            new_list = []
            for image in all_wavelength_dict[self.wavelength]['full_raw_image_data']:
                new_list.append(strided_rescale(image, 8))
                all_wavelength_dict[self.wavelength]['raw_image_data'] = new_list
            image_list = all_wavelength_dict[self.wavelength]['raw_image_data']
            self.image_list = make_HMI_image_list(image_list, -50, 50)
        else:
            self.image_list = all_wavelength_dict[self.wavelength][self.get_image_size()]

    def get_image_size(self):
        if self.size == 'full':
            return 'image_data'
        elif self.size == 'cropped':
            return 'cropped_image_data'

    def draw_RHESSI_limb(self):
        self.limb_state = False
        img = Image.open('real_limb.png')
        sun_radius = all_wavelength_dict[self.wavelength]['sun_radius']
        val = int((sun_radius / all_wavelength_dict[self.wavelength][
            'pixel_arcsec_x']) * 2)  # divide 8 to lower resolution and multiply 2 cuz diamter = radius*2

        img = img.resize((val, val))
        numpydata = asarray(img)
        self.limb_data = np.rot90(numpydata)
        self.limb_x = int((all_wavelength_dict[self.wavelength]['sun_center_x'] - (
                    sun_radius / all_wavelength_dict[self.wavelength]['pixel_arcsec_x'])))
        self.limb_y = int((all_wavelength_dict[self.wavelength]['sun_center_y'] - (
                    sun_radius / all_wavelength_dict[self.wavelength]['pixel_arcsec_y'])))

    def initialize_menu_slider(self):
        self.image_slider_group = QtWidgets.QFrame()
        self.image_slider_group.setGeometry(QtCore.QRect(0, 0, 10, 10))
        self.image_slider_group.setObjectName("RHESSI_group")
        self.image_slider_gridLayout = QtWidgets.QGridLayout(self.image_slider_group)
        self.image_slider_text = QtWidgets.QLabel(self.image_slider_group)
        self.image_slider_text.setObjectName("RHESSI_300")
        self.image_slider_text.setText("Min                                       Max")
        self.image_slider_gridLayout.addWidget(self.image_slider_text, 0, 0, 1, 1)
        self.image_range_slider = QRangeSlider(self.image_slider_group)
        self.image_range_slider.setObjectName("speed_slider")
        self.image_range_slider.setBackgroundStyle('background: black')
        self.image_range_slider.handle.setStyleSheet('background: white')
        if self.wavelength == 'HMI':
            self.image_range_slider.setMin(-1000)
            self.image_range_slider.setMax(1000)
            self.image_range_slider.setRange(-50, 50)
        else:
            self.image_range_slider.setMin(0)
            self.image_range_slider.setMax(self.max_clip_val * 3)
            self.image_range_slider.setRange(0, self.max_clip_val)
        self.image_slider_gridLayout.addWidget(self.image_range_slider, 1, 0, 1, 1)
        self.image_range_slider.startValueChanged.connect(self.image_range_slider_func)
        self.image_range_slider.endValueChanged.connect(self.image_range_slider_func)
        self.time_string_list = all_wavelength_dict[self.wavelength]['time_string_list']

        image_slider_menuItems = [
            (translate("PlotItem", 'Image Limit Slider'), self.image_slider_group)]

        for name, grp in image_slider_menuItems:
            menu = self.view.menu.addMenu(name)
            menu.setStyleSheet("QMenu::item { height: 1000px;margin: 0px;}")
            act = QtWidgets.QWidgetAction(menu)
            act.setDefaultWidget(grp)
            menu.addAction(act)


    def clip_image_list(self, image_list, min_val, max_val):
        new_list = []
        for image in image_list:
            image1 = np.clip(image, min_val, max_val)
            new_list.append(image1)
        return new_list

    def image_range_slider_func(self):
        min_val = self.image_range_slider.start()
        max_val = self.image_range_slider.end()

        self.max_clip_val = max_val

        if self.wavelength == 'HMI' and self.size == 'full':
            self.image_list = make_HMI_image_list(all_wavelength_dict[self.wavelength]['raw_image_data'], min_val,
                                                  max_val)

        elif self.wavelength == 'HMI' and self.size == 'cropped':
            self.image_list = self.clip_image_list(all_wavelength_dict[self.wavelength]['raw_cropped_image_data'],
                                                   min_val,
                                                   max_val)

        elif self.size == 'cropped':
            sdoaia = ct.aia_color_table(self.wavelength * u.angstrom)
            self.image_list = asinh_stretch_images(sdoaia, min_val, max_val,
                                                   all_wavelength_dict[self.wavelength]['raw_cropped_image_data'])
        else:
            sdoaia = ct.aia_color_table(self.wavelength * u.angstrom)
            self.image_list = asinh_stretch_images(sdoaia, min_val, max_val,
                                                   all_wavelength_dict[self.wavelength]['raw_image_data'])

    def change_cropped_aia_map_axis(self):
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
        small_x = changing_solar_movies_dict["small_x"]
        small_y = changing_solar_movies_dict["small_y"]

        if (rounded_small_x - helioprojective_small_x) < 59:
            fourth_x = x_helioprojective_to_pixel(rounded_small_x, 300) - small_x
            x_labels = [(first_x, str(rounded_small_x) + '"'), (second_x, str(rounded_small_x + 100) + '"'),
                        (third_x, str(rounded_small_x + 200) + '"'), (fourth_x, str(rounded_small_x + 300) + '"')]
        else:
            x_labels = [(first_x, str(rounded_small_x) + '"'), (second_x, str(rounded_small_x + 100) + '"'),
                        (third_x, str(rounded_small_x + 200) + '"')]

        if (rounded_small_y - helioprojective_small_y) < 1000:
            fourth_y = y_helioprojective_to_pixel(rounded_small_y, 300) - small_y
            y_labels = [(first_y, str(rounded_small_y) + '"'), (second_y, str(rounded_small_y + 100) + '"'),
                        (third_y, str(rounded_small_y + 200) + '"'), (fourth_y, str(rounded_small_y + 300) + '"')]

        else:
            y_labels = [(first_y, str(rounded_small_y) + '"'), (second_y, str(rounded_small_y + 100) + '"'),
                        (third_y, str(rounded_small_y + 200) + '"')]

        self.axis_left.setTicks([y_labels])
        self.axis_bottom.setTicks([x_labels])

    def change_cropped_RHESSI_map_axis(self, x, y):
        self.view.setLimits(xMin=0, xMax=120, yMin=0, yMax=120)
        square_side_length = 15
        small_y = (512 - y - square_side_length) * 8
        small_x = (x - square_side_length) * 8

        helioprojective_small_x = x_pixel_to_helioprojective(small_x, sun_center_x_list_placeholder[0], pixel_arcsec_x_list_placeholder[0])
        helioprojective_small_y = y_pixel_to_helioprojective(small_y, sun_center_y_list_placeholder[0], pixel_arcsec_y_list_placeholder[0])
        rounded_small_x = roundup(helioprojective_small_x)
        rounded_small_y = rounddown(helioprojective_small_y)

        distance = 200

        first_x = (x_helioprojective_to_pixel(rounded_small_x, 0) - small_x) / 8
        second_x = (x_helioprojective_to_pixel(rounded_small_x, distance) - small_x) / 8
        third_x = (x_helioprojective_to_pixel(rounded_small_x, distance * 2) - small_x) / 8
        first_y = (y_helioprojective_to_pixel(rounded_small_y, 0) - small_y) / 8
        second_y = (y_helioprojective_to_pixel(rounded_small_y, distance) - small_y) / 8
        third_y = (y_helioprojective_to_pixel(rounded_small_y, distance * 2) - small_y) / 8
        fourth_x = (x_helioprojective_to_pixel(rounded_small_x, distance * 3) - small_x) / 8
        fourth_y = (y_helioprojective_to_pixel(rounded_small_y, distance * 3) - small_y) / 8

        if (rounded_small_x - helioprojective_small_x) < 10:
            x_labels = [(first_x, str(rounded_small_x) + "'"), (second_x, str(rounded_small_x + distance) + "'"),
                        (third_x, str(rounded_small_x + distance * 2) + "'"),
                        (fourth_x, str(rounded_small_x + distance * 3) + "'")]
        else:
            x_labels = [(first_x, str(rounded_small_x) + "'"), (second_x, str(rounded_small_x + distance) + "'"),
                        (third_x, str(rounded_small_x + distance * 2) + "'"),
                        (fourth_x, str(rounded_small_x + distance * 3) + "'")]

        if (rounded_small_y - helioprojective_small_y) < 50:
            y_labels = [(first_y, str(rounded_small_y) + "'"), (second_y, str(rounded_small_y + distance) + "'"),
                        (third_y, str(rounded_small_y + distance * 2) + "'"),
                        (fourth_y, str(rounded_small_y + distance * 3) + "'")]

        else:
            y_labels = [(first_y, str(rounded_small_y) + '"'), (second_y, str(rounded_small_y + distance) + '"'),
                        (third_y, str(rounded_small_y + distance * 2) + '"'),
                        (fourth_y, str(rounded_small_y + distance * 3) + '"')]

        self.axis_left.setTicks([y_labels])
        self.axis_bottom.setTicks([x_labels])

    def create_plot(self):
        self.ax2D = self.win.addPlot()

    def draw_limb(self):
        self.limb_state = False
        img = Image.open('real_limb.png')
        sun_radius = all_wavelength_dict[self.wavelength]['sun_radius']
        val = int((sun_radius / all_wavelength_dict[self.wavelength][
            'pixel_arcsec_x']) / 8 * 2)  # divide 8 to lower resolution and multiply 2 cuz diamter = radius*2

        img = img.resize((val, val))
        numpydata = asarray(img)
        self.limb_data = np.rot90(numpydata)
        self.limb_x = int((all_wavelength_dict[self.wavelength]['sun_center_x'] - (
                    sun_radius / all_wavelength_dict[self.wavelength]['pixel_arcsec_x'])) / 8)
        self.limb_y = int((all_wavelength_dict[self.wavelength]['sun_center_y'] - (
                    sun_radius / all_wavelength_dict[self.wavelength]['pixel_arcsec_y'])) / 8)

    def putLabel(self, label, axes, x=110, y=60):
        label.setPos(QtCore.QPointF(x, y))
        label.setFont(QFont("Times New Roman", 22))
        label.setFont(QFont("Times New Roman", 22))
        axes.addItem(label)

    def put_wavelength_label(self, label, axes, x=30, y=485):
        label.setPos(QtCore.QPointF(x, y))
        label.setFont(QFont("Times New Roman", 15))
        axes.addItem(label)
        label.setText('AIA ' + str(self.wavelength) + ' Å')

    def update_image_list(self):
        state_to_list_dict = {'Normal': self.get_image_size(), 'Running Difference': 'run_diff_list',
                              'Base Difference': 'base_diff_list',
                              'Running Ratio': 'run_ratio_list', 'Base Ratio': 'base_ratio_list'}
        self.image_list = all_wavelength_dict[self.wavelength][state_to_list_dict[self.state]]

    def aia_map_change_wavelength_func(self):
        a = self.window.sender().text()
        print(a)
        if a == 'HMI' or a == 'RHESSI':
            self.wavelength = a
        else:
            self.wavelength = int(a)
        self.wavelength_label.setText('AIA ' + str(self.wavelength) + ' Å')
        value_list = list(self.without_keys(self.wavelength_menu_dict, str(self.wavelength)).values())
        for value in value_list:
            value.setChecked(False)
        self.time_string_list = all_wavelength_dict[self.wavelength]['time_string_list']
        if self.wavelength == 'RHESSI':
            self.change_cropped_RHESSI_map_axis(RHESSI_center_x, RHESSI_center_y)
        else:
            self.view.setLimits(xMin=0, xMax=600, yMin=0, yMax=600)
            self.change_cropped_aia_map_axis()
        self.update_image_list()

    def aia_map_change_image_func(self):
        a = self.window.sender().text()
        self.state = str(a)
        value_list = list(self.without_keys(self.ratio_menu_dict, self.state).values())
        for value in value_list:
            value.setChecked(False)
        self.update_image_list()

    def update_aia_image(self, i):
        try:
            self.aia_image.setImage(self.image_list[i])

        except IndexError:
            self.aia_image.setImage(self.image_list[0])

    def hide_limb(self):
        self.limb_state = False
        self.ax2D.removeItem(self.limb_img)

    def show_limb_func(self):
        self.limb_state = True
        self.limb_img = pg.ImageItem(view=self.plot)
        self.limb_img.setImage(self.limb_data)
        self.limb_img.setPos(QtCore.QPointF(self.limb_x, self.limb_y))
        self.ax2D.addItem(self.limb_img)

    def put_limb_image(self):
        print(self.limb_state)
        if self.limb_state == True:
            self.hide_limb()
        else:
            self.show_limb_func()

    def display_func(self):
        self.ax2D.setAspectLocked(True)
        self.axis_bottom = self.ax2D.getAxis('bottom')
        self.axis_left = self.ax2D.getAxis('left')
        self.view = self.ax2D.getViewBox()
        self.axis_left.setZValue(1)
        self.axis_bottom.setZValue(1)
        self.axis_bottom.setTickSpacing(100, 100)
        self.axis_left.setTickSpacing(100, 100)
        self.view.setBorder()
        self.aia_image = pg.ImageItem(border='w', view=self.plot)
        self.ax2D.addItem(self.aia_image)
        self.putLabel(self.time_label, self.ax2D)
        self.put_wavelength_label(self.wavelength_label, self.ax2D)
        self.show_limb = self.view.menu.addAction('Show Limb')
        self.show_limb.setCheckable(True)
        self.show_limb.setChecked(True)
        self.show_limb.triggered.connect(self.put_limb_image)

        # self.put_limb_image(plot, self.ax2D, x=self.limb_x, y=self.limb_y)

    def without_keys(self, d, keys):
        return {x: d[x] for x in d if x not in keys}

    def setup_wavelength_button(self, button):
        button.setCheckable(True)
        button.triggered.connect(self.aia_map_change_wavelength_func)

    def initialize_wavelength_menu(self):
        wavelength_menu = self.view.menu.addMenu('Change Wavelength')
        self.change_to_94 = wavelength_menu.addAction('94')
        self.change_to_171 = wavelength_menu.addAction('171')
        self.change_to_211 = wavelength_menu.addAction('211')
        self.change_to_304 = wavelength_menu.addAction('304')

        self.change_to_1600 = wavelength_menu.addAction('1600')
        self.change_to_RHESSI = wavelength_menu.addAction('RHESSI')
        self.wavelength_menu_dict = {'94': self.change_to_94, '171': self.change_to_171,'211': self.change_to_211,
                                    '304': self.change_to_304,
                                     '1600': self.change_to_1600, 'RHESSI': self.change_to_RHESSI}
        for button in self.wavelength_menu_dict.values():
            self.setup_wavelength_button(button)

    def setup_ratio_button(self, button):
        button.setCheckable(True)
        button.triggered.connect(self.aia_map_change_image_func)

    def initialize_ratio_menu(self):
        difference_menu = self.view.menu.addMenu('Change Display')
        self.normal = difference_menu.addAction('Normal')
        self.running_diff = difference_menu.addAction('Running Difference')
        self.base_diff = difference_menu.addAction('Base Difference')
        self.running_ratio = difference_menu.addAction('Running Ratio')
        self.base_ratio = difference_menu.addAction('Base Ratio')

        self.ratio_menu_dict = {'Normal': self.normal, 'Running Difference': self.running_diff,
                                'Base Difference': self.base_diff,
                                'Running Ratio': self.running_ratio, 'Base Ratio': self.base_ratio}

        for button in self.ratio_menu_dict.values():
            self.setup_ratio_button(button)

        self.normal.setChecked(True)

    def show_center_mark(self, i):
        try:
            lineitem = self.show_center_list[0]
            self.ax2D.removeItem(lineitem)
            lineitem2 = self.show_center_list[1]
            self.ax2D.removeItem(lineitem2)
        except (AttributeError, IndexError) as err:
            pass

        self.show_center_list = []
        x, y = brightest_pixel_dict[i]['cropped_pixel']
        self.line = LineItem(QtCore.QLineF(x - 20, y - 20, x + 20, y + 20))
        self.line.setPen(QtCore.Qt.red)
        self.ax2D.addItem(self.line)
        self.show_center_list.append(self.line)

        self.line2 = LineItem(QtCore.QLineF(x - 20, y + 20, x + 20, y - 20))
        pen = QPen()  # creates a default pen
        pen.setWidth(3)
        pen.setBrush(Qt.red)
        self.line.setPen(pen)
        self.line2.setPen(pen)
        self.ax2D.addItem(self.line2)
        self.show_center_list.append(self.line2)

    def show_center_func(self):
        if self.show_center == True:
            self.show_center = False
            try:
                lineitem = self.show_center_list[0]
                self.ax2D.removeItem(lineitem)
                lineitem2 = self.show_center_list[1]
                self.ax2D.removeItem(lineitem2)
            except (AttributeError, IndexError) as err:
                pass
        else:
            self.show_center = True

    def initialize_show_center_button(self):
        self.show_center = False
        self.show_center_button = self.view.menu.addAction('Show Center')
        self.show_center_button.setCheckable(True)
        self.show_center_button.triggered.connect(self.show_center_func)

    def update_axis_dict(self, ix, iy):
        square_side_length = 300
        small_x = ix - square_side_length
        small_y = iy - square_side_length
        big_x = ix + square_side_length
        big_y = iy + square_side_length

        helioprojective_small_x = x_pixel_to_helioprojective(small_x, sun_center_x_list_placeholder[0], pixel_arcsec_x_list_placeholder[0])
        helioprojective_small_y = y_pixel_to_helioprojective(small_y, sun_center_y_list_placeholder[0], pixel_arcsec_y_list_placeholder[0])
        rounded_small_x = roundup(helioprojective_small_x)
        rounded_small_y = rounddown(helioprojective_small_y)

        first_x = x_helioprojective_to_pixel(rounded_small_x, 0) - small_x
        second_x = x_helioprojective_to_pixel(rounded_small_x, 100) - small_x
        third_x = x_helioprojective_to_pixel(rounded_small_x, 200) - small_x
        first_y = y_helioprojective_to_pixel(rounded_small_y, 0) - small_y
        second_y = y_helioprojective_to_pixel(rounded_small_y, 100) - small_y
        third_y = y_helioprojective_to_pixel(rounded_small_y, 200) - small_y

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
        user_projective_x = x_pixel_to_helioprojective(ix, sun_center_x_list_placeholder[0], pixel_arcsec_x_list_placeholder[0])
        user_projective_y = y_pixel_to_helioprojective(iy, sun_center_y_list_placeholder[0], pixel_arcsec_y_list_placeholder[0])
        changing_solar_movies_dict[
            "user_selected_center_helioprojective"] = user_projective_x, user_projective_y

        changing_solar_movies_dict["user_selected_center_heliographical"] = helioprojective_to_heliographic(
            user_projective_x, user_projective_y, aia_map_object)

    def set_axis(self):
        self.view.setLimits(xMin=0, xMax=512, yMin=0, yMax=512)
        x_labels = [((self.sun_center_x - (1000 / self.pixel_arcsec_x)) / 8, "-1000'"),
                    ((self.sun_center_x - (500 / self.pixel_arcsec_x)) / 8, "-500'"), (self.sun_center_x / 8, "0'"),
                    ((self.sun_center_x + (500 / self.pixel_arcsec_x)) / 8, "500'"),
                    ((self.sun_center_x + (1000 / self.pixel_arcsec_x)) / 8, "1000'")]

        y_labels = [((self.sun_center_y - (1000 / self.pixel_arcsec_y)) / 8, "-1000'"),
                    ((self.sun_center_y - (500 / self.pixel_arcsec_y)) / 8, "-500'"), (self.sun_center_y / 8, "0'"),
                    ((self.sun_center_y + (500 / self.pixel_arcsec_y)) / 8, "500'"),
                    ((self.sun_center_y + (1000 / self.pixel_arcsec_y)) / 8, "1000'")]

        self.axis_bottom.setTicks([x_labels])
        self.axis_left.setTicks([y_labels])

    def draw_rectangle(self, x, y):
        # rect_item.setFlag(QtWidgets.QGraphicsItem.ItemIsMovable, True)
        try:
            self.ax2D.removeItem(self.rect)
        except Exception:
            pass

        self.rect = RectItem(QtCore.QRectF(x / 8 - 37.5, y / 8 - 37.5, 75, 75))
        self.rect.setPen(QtCore.Qt.red)
        self.ax2D.addItem(self.rect)

    def change_image_clip_values(self, clip_min_value, clip_max_value):
        sdoaia = ct.aia_color_table(self.wavelength * u.angstrom)
        image_dict = all_wavelength_dict[self.wavelength]

        run_diff_list, base_diff_list, run_ratio_list, base_ratio_list = asinh_stretch_images(sdoaia, clip_min_value,
                                                                                              clip_max_value,
                                                                                              image_dict[
                                                                                                  'run_diff_list'],
                                                                                              image_dict[
                                                                                                  'base_diff_list'],
                                                                                              image_dict[
                                                                                                  'run_ratio_list'],
                                                                                              image_dict[
                                                                                                  'base_ratio_list'])
        all_wavelength_dict[self.wavelength]['run_diff_list'] = run_diff_list
        all_wavelength_dict[self.wavelength]['base_diff_list'] = base_diff_list
        all_wavelength_dict[self.wavelength]['run_ratio_list'] = run_ratio_list
        all_wavelength_dict[self.wavelength]['base_ratio_list'] = base_ratio_list


class Full_RHESSI_Widget(AIA_Image_Widgets):
    def __init__(self, win, window, wavelength, size):
        super().__init__(win, window, wavelength, size)
        self.create_plot()
        self.display_func()
        self.set_RHESSI_axis()
        self.draw_RHESSI_limb()
        self.draw_image()
        self.show_limb.setChecked(False)
        self.time_label.setText(all_wavelength_dict['RHESSI']['time_string_list'][0])
        self.wavelength_label.setText("RHESSI")
        self.putLabel(self.time_label, self.ax2D, x=30, y=60)
        self.put_wavelength_label(self.wavelength_label, self.ax2D)
        self.draw_rectangle(RHESSI_center_x * 8, RHESSI_center_y * 8)

    def draw_image(self):
        self.aia_image.setImage(all_wavelength_dict[self.wavelength]['image_data'])

    def set_RHESSI_axis(self):
        self.view.setLimits(xMin=0, xMax=512, yMin=0, yMax=512)
        x_labels = [((self.sun_center_x - (1000 / self.pixel_arcsec_x)), "-1000'"),
                    ((self.sun_center_x - (500 / self.pixel_arcsec_x)), "-500'"), (self.sun_center_x, "0'"),
                    ((self.sun_center_x + (500 / self.pixel_arcsec_x)), "500'"),
                    ((self.sun_center_x + (1000 / self.pixel_arcsec_x)), "1000'")]

        y_labels = [((self.sun_center_y - (1000 / self.pixel_arcsec_y)), "-1000'"),
                    ((self.sun_center_y - (500 / self.pixel_arcsec_y)), "-500'"), (self.sun_center_y, "0'"),
                    ((self.sun_center_y + (500 / self.pixel_arcsec_y)), "500'"),
                    ((self.sun_center_y + (1000 / self.pixel_arcsec_y)), "1000'")]

        self.axis_bottom.setTicks([x_labels])
        self.axis_left.setTicks([y_labels])


class Full_Helioprojective_HMI_Widget(AIA_Image_Widgets):
    def __init__(self, win, window, wavelength, size):
        super().__init__(win, window, wavelength, size)
        self.create_plot()
        self.display_func()
        self.draw_limb()
        self.set_axis()
        self.initialize_menu_slider()
        self.draw_rectangle(peak_x, peak_y)
        self.wavelength_label.setPos(QtCore.QPointF(30, 485))
        self.wavelength_label.setFont(QFont("Times New Roman", 15))
        self.ax2D.addItem(self.wavelength_label)
        self.wavelength_label.setText('HMI Magnetogram')
        self.show_limb.setChecked(False)


class Other_AIA_Widgets(AIA_Image_Widgets):
    def __init__(self, win, window, wavelength, size):
        super().__init__(win, window, wavelength, size)
        self.create_plot()
        self.display_func()
        self.initialize_menu_slider()
        self.initialize_wavelength_menu()
        self.initialize_ratio_menu()
        self.set_axis()
        self.draw_limb()
        self.put_limb_image()
        self.draw_rectangle(peak_x, peak_y)


class Full_AIA_Map_Widget(pg.PlotWidget):
    def __init__(self, wavelength, window, size1, ui_active_region_number, ui_hale_class_value, cropped_map_object,
                 full_hmi_map_object, cropped_1600_map_object, cropped_hmi_map_object, full_1600_map_object,
                 full_helioprojective_hmi_map_object, graph_list, normalized_graph_list):
        pg.PlotWidget.__init__(self)
        pg.setConfigOptions(antialias=True)
        self.window = window
        self.time_label = pg.TextItem('helo', **{'color': '#ffffff'})
        self.wavelength_label = pg.TextItem('helo', **{'color': '#ffffff'})
        self.wavelength = wavelength
        self.time_string_list = all_wavelength_dict[self.wavelength]['time_string_list']
        # self.max_clip_val = all_wavelength_dict[self.wavelength]['max_clip_val']
        self.sun_center_x = sun_center_x_list_placeholder[0]
        self.sun_center_y = sun_center_y_list_placeholder[0]
        self.pixel_arcsec_x = pixel_arcsec_x_list_placeholder[0]
        self.pixel_arcsec_y = pixel_arcsec_y_list_placeholder[0]
        self.state = 'Normal'
        self.size1 = size1
        self.max_clip_val = int(all_wavelength_dict[self.wavelength]['max_clip_val'])
        self.image_list = all_wavelength_dict[self.wavelength][self.get_image_size()]

        self.active_region_number = ui_active_region_number
        self.hale_class_value = ui_hale_class_value

        self.view = self.getViewBox()
        self.display_func()
        self.initialize_menu_slider()
        self.initialize_wavelength_menu()
        self.initialize_ratio_menu()
        self.set_axis()
        self.draw_rectangle(peak_x, peak_y)
        self.draw_limb()
        self.put_limb_image()
        self.AR_number_labels(self)
        self.cropped_hmi_object = cropped_hmi_map_object
        self.cropped_map_object = cropped_map_object
        self.full_hmi_image = full_hmi_map_object
        self.cropped_1600_object = cropped_1600_map_object
        self.full_helioprojective_HMI_object = full_helioprojective_hmi_map_object
        self.full_1600_map_object = full_1600_map_object
        self.cropped_widget_list = [self.cropped_map_object, self.cropped_1600_object, self.cropped_hmi_object]
        self.graph_list = graph_list
        self.normalized_graph_list = normalized_graph_list

    def update_axis_dict(self, ix, iy):
        square_side_length = 300
        small_x = ix - square_side_length
        small_y = iy - square_side_length
        big_x = ix + square_side_length
        big_y = iy + square_side_length

        helioprojective_small_x = x_pixel_to_helioprojective(small_x, sun_center_x_list_placeholder[0], pixel_arcsec_x_list_placeholder[0])
        helioprojective_small_y = y_pixel_to_helioprojective(small_y, sun_center_y_list_placeholder[0], pixel_arcsec_y_list_placeholder[0])
        rounded_small_x = roundup(helioprojective_small_x)
        rounded_small_y = rounddown(helioprojective_small_y)

        first_x = x_helioprojective_to_pixel(rounded_small_x, 0) - small_x
        second_x = x_helioprojective_to_pixel(rounded_small_x, 100) - small_x
        third_x = x_helioprojective_to_pixel(rounded_small_x, 200) - small_x
        first_y = y_helioprojective_to_pixel(rounded_small_y, 0) - small_y
        second_y = y_helioprojective_to_pixel(rounded_small_y, 100) - small_y
        third_y = y_helioprojective_to_pixel(rounded_small_y, 200) - small_y

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
        user_projective_x = x_pixel_to_helioprojective(ix,  sun_center_x_list_placeholder[0], pixel_arcsec_x_list_placeholder[0])
        user_projective_y = y_pixel_to_helioprojective(iy, sun_center_y_list_placeholder[0], pixel_arcsec_y_list_placeholder[0])
        changing_solar_movies_dict[
            "user_selected_center_helioprojective"] = user_projective_x, user_projective_y

        changing_solar_movies_dict["user_selected_center_heliographical"] = helioprojective_to_heliographic(
            user_projective_x, user_projective_y, aia_map_object)

    def display_func(self):
        self.setAspectLocked(True)
        self.axis_bottom = self.getAxis('bottom')
        self.axis_left = self.getAxis('left')
        self.view = self.getViewBox()
        self.axis_left.setZValue(1)
        self.axis_bottom.setZValue(1)
        self.axis_bottom.setTickSpacing(100, 100)
        self.axis_left.setTickSpacing(100, 100)
        self.view.setBorder()
        self.aia_image = pg.ImageItem(border='w')
        self.addItem(self.aia_image)
        self.putLabel(self.time_label, self)
        self.put_wavelength_label(self.wavelength_label, self)

        show_limb = self.view.menu.addAction('Show Limb')

        show_limb.setCheckable(True)
        show_limb.setChecked(True)
        show_limb.triggered.connect(self.put_limb_image)

    def putLabel(self, label, axes, x=110, y=60):
        label.setPos(QtCore.QPointF(x, y))
        label.setFont(QFont("Times New Roman", 22))
        label.setFont(QFont("Times New Roman", 22))
        axes.addItem(label)

    def put_wavelength_label(self, label, axes, x=30, y=485):
        label.setPos(QtCore.QPointF(x, y))
        label.setFont(QFont("Times New Roman", 15))
        axes.addItem(label)
        label.setText('AIA ' + str(self.wavelength) + ' Å')

    def get_image_size(self):
        if self.size1 == 'full':
            return 'image_data'
        elif self.size1 == 'cropped':
            return 'cropped_image_data'

    def initialize_menu_slider(self):
        self.image_slider_group = QtWidgets.QFrame()
        self.image_slider_group.setGeometry(QtCore.QRect(0, 0, 10, 10))
        self.image_slider_group.setObjectName("RHESSI_group")
        self.image_slider_gridLayout = QtWidgets.QGridLayout(self.image_slider_group)
        self.image_slider_text = QtWidgets.QLabel(self.image_slider_group)
        self.image_slider_text.setObjectName("RHESSI_300")
        self.image_slider_text.setText("Min                                       Max")
        self.image_slider_gridLayout.addWidget(self.image_slider_text, 0, 0, 1, 1)
        self.image_range_slider = QRangeSlider(self.image_slider_group)
        self.image_range_slider.setObjectName("speed_slider")
        self.image_range_slider.setBackgroundStyle('background: black')
        self.image_range_slider.handle.setStyleSheet('background: white')
        self.image_range_slider.setMin(0)
        self.image_range_slider.setMax(self.max_clip_val*3)
        self.image_range_slider.setRange(0, int(self.max_clip_val))
        self.image_slider_gridLayout.addWidget(self.image_range_slider, 1, 0, 1, 1)
        self.image_range_slider.startValueChanged.connect(self.image_range_slider_func)
        self.image_range_slider.endValueChanged.connect(self.image_range_slider_func)
        self.time_string_list = all_wavelength_dict[self.wavelength]['time_string_list']

        image_slider_menuItems = [
            (translate("PlotItem", 'Image Limit Slider'), self.image_slider_group)]

        for name, grp in image_slider_menuItems:
            menu = self.view.menu.addMenu(name)
            menu.setStyleSheet("QMenu::item { height: 1000px;margin: 0px;}")
            act = QtWidgets.QWidgetAction(menu)
            act.setDefaultWidget(grp)
            menu.addAction(act)

    def image_range_slider_func(self):
        min_val = self.image_range_slider.start()
        max_val = self.image_range_slider.end()

        if self.wavelength == 'HMI' and self.size == 'full':
            self.image_list = make_HMI_image_list(all_wavelength_dict[self.wavelength]['raw_image_data'], min_val,
                                                  max_val)
        elif self.size == 'cropped':
            self.image_list = self.clip_image_list(all_wavelength_dict[self.wavelength]['cropped_raw_image_data'],
                                                   min_val, max_val)
        else:
            sdoaia = ct.aia_color_table(self.wavelength * u.angstrom)
            self.image_list = asinh_stretch_images(sdoaia, min_val, max_val,
                                                   all_wavelength_dict[self.wavelength]['raw_image_data'])

    def without_keys(self, d, keys):
        return {x: d[x] for x in d if x not in keys}

    def setup_wavelength_button(self, button):
        button.setCheckable(True)
        button.triggered.connect(self.aia_map_change_wavelength_func)

    def initialize_wavelength_menu(self):
        wavelength_menu = self.view.menu.addMenu('Change Wavelength')
        self.change_to_94 = wavelength_menu.addAction('94')
        self.change_to_211 = wavelength_menu.addAction('211')
        self.change_to_171 = wavelength_menu.addAction('171')
        self.change_to_304 = wavelength_menu.addAction('304')
        self.change_to_1600 = wavelength_menu.addAction('1600')

        self.wavelength_menu_dict = {'94': self.change_to_94, '171': self.change_to_171,'211': self.change_to_211,'304': self.change_to_304}
        for button in self.wavelength_menu_dict.values():
            self.setup_wavelength_button(button)

    def setup_ratio_button(self, button):
        button.setCheckable(True)
        button.triggered.connect(self.aia_map_change_image_func)

    def initialize_ratio_menu(self):
        difference_menu = self.view.menu.addMenu('Change Display')
        self.normal = difference_menu.addAction('Normal')
        self.running_diff = difference_menu.addAction('Running Difference')
        self.base_diff = difference_menu.addAction('Base Difference')
        self.running_ratio = difference_menu.addAction('Running Ratio')
        self.base_ratio = difference_menu.addAction('Base Ratio')

        self.ratio_menu_dict = {'Normal': self.normal, 'Running Difference': self.running_diff,
                                'Base Difference': self.base_diff,
                                'Running Ratio': self.running_ratio, 'Base Ratio': self.base_ratio}

        for button in self.ratio_menu_dict.values():
            self.setup_ratio_button(button)

        self.normal.setChecked(True)

    def show_center_mark(self, i):
        try:
            lineitem = self.show_center_list[0]
            self.removeItem(lineitem)
            lineitem2 = self.show_center_list[1]
            self.removeItem(lineitem2)
        except (AttributeError, IndexError) as err:
            pass

        self.show_center_list = []
        x, y = brightest_pixel_dict[i]['cropped_pixel']
        self.line = LineItem(QtCore.QLineF(x - 20, y - 20, x + 20, y + 20))
        self.line.setPen(QtCore.Qt.red)
        self.addItem(self.line)
        self.show_center_list.append(self.line)

        self.line2 = LineItem(QtCore.QLineF(x - 20, y + 20, x + 20, y - 20))
        pen = QPen()
        pen.setWidth(3)
        pen.setBrush(Qt.red)
        self.line.setPen(pen)
        self.line2.setPen(pen)
        self.addItem(self.line2)
        self.show_center_list.append(self.line2)

    def show_center_func(self):
        if self.show_center == True:
            self.show_center = False
            try:
                lineitem = self.show_center_list[0]
                self.removeItem(lineitem)
                lineitem2 = self.show_center_list[1]
                self.removeItem(lineitem2)
            except (AttributeError, IndexError) as err:
                pass
        else:
            self.show_center = True

    def initialize_show_center_button(self):
        self.show_center = False
        self.show_center_button = self.view.menu.addAction('Show Center')
        self.show_center_button.triggered.connect(self.show_center_func)

    def set_axis(self):
        self.view.setLimits(xMin=0, xMax=512, yMin=0, yMax=512)
        x_labels = [((self.sun_center_x - (1000 / self.pixel_arcsec_x)) / 8, "-1000'"),
                    ((self.sun_center_x - (500 / self.pixel_arcsec_x)) / 8, "-500'"), (self.sun_center_x / 8, "0'"),
                    ((self.sun_center_x + (500 / self.pixel_arcsec_x)) / 8, "500'"),
                    ((self.sun_center_x + (1000 / self.pixel_arcsec_x)) / 8, "1000'")]

        y_labels = [((self.sun_center_y - (1000 / self.pixel_arcsec_y)) / 8, "-1000'"),
                    ((self.sun_center_y - (500 / self.pixel_arcsec_y)) / 8, "-500'"), (self.sun_center_y / 8, "0'"),
                    ((self.sun_center_y + (500 / self.pixel_arcsec_y)) / 8, "500'"),
                    ((self.sun_center_y + (1000 / self.pixel_arcsec_y)) / 8, "1000'")]

        self.axis_bottom.setTicks([x_labels])
        self.axis_left.setTicks([y_labels])

    def draw_rectangle(self, x, y):
        # rect_item.setFlag(QtWidgets.QGraphicsItem.ItemIsMovable, True)
        try:
            self.removeItem(self.rect)
        except Exception:
            pass

        self.rect = RectItem(QtCore.QRectF(x / 8 - 37.5, (y) / 8 - 37.5, 75, 75))
        self.rect.setPen(QtCore.Qt.red)
        self.addItem(self.rect)

    def change_image_clip_values(self, clip_min_value, clip_max_value):
        sdoaia = ct.aia_color_table(self.wavelength * u.angstrom)
        image_dict = all_wavelength_dict[self.wavelength]

        run_diff_list, base_diff_list, run_ratio_list, base_ratio_list = asinh_stretch_images(sdoaia, clip_min_value,
                                                                                              clip_max_value,
                                                                                              image_dict[
                                                                                                  'run_diff_list'],
                                                                                              image_dict[
                                                                                                  'base_diff_list'],
                                                                                              image_dict[
                                                                                                  'run_ratio_list'],
                                                                                              image_dict[
                                                                                                  'base_ratio_list'])
        all_wavelength_dict[self.wavelength]['run_diff_list'] = run_diff_list
        all_wavelength_dict[self.wavelength]['base_diff_list'] = base_diff_list
        all_wavelength_dict[self.wavelength]['run_ratio_list'] = run_ratio_list
        all_wavelength_dict[self.wavelength]['base_ratio_list'] = base_ratio_list

    def draw_limb(self):
        self.limb_state = False
        img = Image.open('real_limb.png')
        sun_radius = all_wavelength_dict[self.wavelength]['sun_radius']
        val = int((sun_radius / all_wavelength_dict[self.wavelength][
            'pixel_arcsec_x']) / 8 * 2)  # divide 8 to lower resolution and multiply 2 cuz diamter = radius*2

        img = img.resize((val, val))
        numpydata = asarray(img)
        self.limb_data = np.rot90(numpydata)
        self.limb_x = int((all_wavelength_dict[self.wavelength]['sun_center_x'] - (
                sun_radius / all_wavelength_dict[self.wavelength]['pixel_arcsec_x'])) / 8)
        self.limb_y = int((all_wavelength_dict[self.wavelength]['sun_center_y'] - (
                sun_radius / all_wavelength_dict[self.wavelength]['pixel_arcsec_y'])) / 8)

    def put_limb_image(self):
        if self.limb_state == True:
            self.limb_state = False
            self.removeItem(self.limb_img)
        else:
            self.limb_state = True
            self.limb_img = pg.ImageItem()
            self.limb_img.setImage(self.limb_data)
            self.limb_img.setPos(QtCore.QPointF(self.limb_x, self.limb_y))
            self.addItem(self.limb_img)

    def AR_put_label(self, label, axes, x, y):
        label.setPos(QtCore.QPointF(x - 30, y))
        label.setFont(QFont("Times New Roman", 10))
        axes.addItem(label)

    def update_image_list(self):
        state_to_list_dict = {'Normal': self.get_image_size(), 'Running Difference': 'run_diff_list',
                              'Base Difference': 'base_diff_list',
                              'Running Ratio': 'run_ratio_list', 'Base Ratio': 'base_ratio_list'}
        self.image_list = all_wavelength_dict[self.wavelength][state_to_list_dict[self.state]]

    def aia_map_change_wavelength_func(self):
        a = self.window.sender().text()
        print(a)
        self.wavelength = int(a)
        self.wavelength_label.setText('AIA ' + str(self.wavelength) + ' Å')
        value_list = list(self.without_keys(self.wavelength_menu_dict, str(self.wavelength)).values())
        for value in value_list:
            value.setChecked(False)
        self.time_string_list = all_wavelength_dict[self.wavelength]['time_string_list']
        self.update_image_list()

    def aia_map_change_image_func(self):
        a = self.window.sender().text()
        print(a)
        self.state = str(a)
        value_list = list(self.without_keys(self.ratio_menu_dict, self.state).values())
        for value in value_list:
            value.setChecked(False)
        self.update_image_list()

    def update_aia_image(self, i):
        try:
            self.aia_image.setImage(self.image_list[i])
        except IndexError:
            self.aia_image.setImage(self.image_list[0])

    def AR_number_labels(self, axes):
        x = peak_x / 8
        y = peak_y / 8
        real_difference = 9999999
        real_AR_number = 0
        hale_class = 'beta'
        for region2 in AR_dict:
            region1 = AR_dict[region2]
            time_label = pg.TextItem('helo', **{'color': '#ffffff'})
            time_label.setText(str(region2))
            region_x = int(region1[0])
            region_y = int(region1[1])

            difference = abs(region_x - x) + abs(region_y - y)
            if real_difference > difference:
                real_difference = difference
                real_AR_number = region1[2]
                hale_class = region1[3]

            self.AR_put_label(time_label, axes, region_x, region_y)

        self.active_region_number.setText(str(real_AR_number))
        self.hale_class_value.setText(str(hale_class))

    def mouseDoubleClickEvent(self, event):
        if event.button() == QtCore.Qt.RightButton:
            pass
        else:
            scene_pos = self.mapToScene(event.pos())
            img_pos = self.aia_image.mapFromScene(scene_pos)
            x_click = img_pos.x()
            y_click = img_pos.y()
            ix = x_click * 8
            iy = y_click * 8
            changing_solar_movies_dict["center"] = [ix, iy]

            user_cropped_map_intensity_list = []

            try:
                self.removeItem(self.user_rect)

            except:
                pass
            self.user_rect = RectItem(QtCore.QRectF(x_click - 37.5, y_click - 37.5, 75, 75))

            self.user_rect.setPen(QtCore.Qt.white)
            self.addItem(self.user_rect)

            self.full_helioprojective_HMI_object.draw_rectangle(ix, iy)
            self.full_1600_map_object.draw_rectangle(ix, iy)

            bruh = 300


            x_click = int(ix)
            y_click = int(iy)

            for wavelength1 in wavelength_list:
                colormap = ct.aia_color_table(wavelength1 * u.angstrom)
                cropped_image_list = []
                intensity_list = []
                for image in all_wavelength_dict[wavelength1]['full_raw_image_data']:
                    cropped_image_data = image[x_click - bruh:x_click + bruh, y_click - bruh:y_click + bruh]
                    cropped_image_list.append(cropped_image_data)

                    intensity_value = np.divide(cropped_image_data.sum(), all_wavelength_dict[wavelength1]['exposure_list'][i])
                    intensity_list.append(intensity_value)
                    normalized_intensity_list = normalized_plot(intensity_list)
                    all_wavelength_dict[wavelength1]["intensity_list"] = intensity_list
                    all_wavelength_dict[wavelength1]["normalized_intensity_list"] = normalized_plot(intensity_list)
                    all_wavelength_dict[wavelength1]['raw_cropped_image_data'] = cropped_image_list
                    self.graph_list[wavelength1].intensity_list = intensity_list
                    self.graph_list[wavelength1].update_data()
                    self.normalized_graph_list[wavelength1].intensity_list = normalized_intensity_list
                    self.normalized_graph_list[wavelength1].update_data()
                    cropped_image_list = asinh_stretch_images(colormap, 0,
                                                              all_wavelength_dict[wavelength1]['max_clip_val'],
                                                              cropped_image_list)
                    all_wavelength_dict[wavelength1]['cropped_image_data'] = cropped_image_list


            self.update_axis_dict(ix, iy)

            self.cropped_map_object.change_cropped_aia_map_axis()
            self.cropped_map_object.image_list =  all_wavelength_dict[self.wavelength]['cropped_image_data']
            self.cropped_hmi_object.create_contours_list(ix, iy)
            self.cropped_1600_object.create_contours_list(ix, iy)

            for widget in self.cropped_widget_list:
                widget.x = int(ix)
                widget.y = int(iy)
                widget.draw_cropped_limb()
                if widget.limb_state == True:
                    widget.hide_limb()
                    widget.show_limb_func()

                if widget.hmi_center_hide_or_show:
                    widget.hide_hmi_center()
                    widget.show_hmi_center()

            brightest_pixel_dict['box_position_pixel'] = [ix, iy]
            helioprojective_x = x_pixel_to_helioprojective(ix, sun_center_x_list_placeholder[0],
                                                           pixel_arcsec_x_list_placeholder[0])
            helioprojective_y = y_pixel_to_helioprojective(iy, sun_center_y_list_placeholder[0],
                                                           pixel_arcsec_y_list_placeholder[0])
            brightest_pixel_dict['box_position_helioprojective'] = [helioprojective_x, helioprojective_y]
            heliographic_coords = helioprojective_to_heliographic(helioprojective_x, helioprojective_y,
                                                                  aia_map_object)
            brightest_pixel_dict['box_position_heliographical'] = heliographic_coords



class Cropped_Maps(AIA_Image_Widgets):
    def __init__(self, win, window, wavelength, size, full_HMI_image):
        super().__init__(win, window, wavelength, size)
        self.create_plot()
        self.display_func()
        self.initialize_menu_slider()
        self.initialize_cropped_image_axis()
        self.change_cropped_aia_map_axis()
        self.initialize_wavelength_menu()
        self.initialize_contour_menu()
        self.initialize_hmi_center_menu()
        self.initialize_polarities_menu()
        self.put_wavelength_label(self.wavelength_label, self.ax2D, x=35, y=568)
        self.contour_display_filter = []
        self.polarity_display_filter = []
        self.full_HMI_image = full_HMI_image
        self.x = peak_x
        self.y = peak_y
        self.initialize_cropped_limb()
        self.draw_cropped_limb()
        self.show_limb.setChecked(False)
        self.limb_state = False
        self.initialize_arrow_menu()

    def initialize_arrow_menu(self):
        CME_angle = flare_metadata['Angle']
        if not CME_angle == 'Halo':
            self.show_arrow_menu = self.view.menu.addAction('Show Arrow')
            self.show_arrow_menu.setCheckable(True)
            self.show_arrow_menu.setChecked(True)
            self.show_arrow_menu.triggered.connect(self.show_arrow_func)
            self.show_arrow_bool = True
            self.show_arrow()

    def show_arrow_func(self):
        if self.show_arrow_bool == True:
            self.hide_arrow()
            self.show_arrow_bool = False
        else:
            self.show_arrow()
            self.show_arrow_bool = True

    def show_arrow(self):
        CME_angle = int(flare_metadata['Angle'])
        if CME_angle == 'Halo':
            pass
        else:
            self.CME_arrow = MyArrowItem(angle=90 - CME_angle, tipAngle=60, headLen=20, tailLen=20, tailWidth=10,
                                         pen={'color': 'w', 'width': 3})
            self.CME_arrow.setPos(512 / 2, 512 / 2)
            self.ax2D.addItem(self.CME_arrow)

    def hide_arrow(self):
        self.ax2D.removeItem(self.CME_arrow)

    def initialize_cropped_limb(self, full_image_val=4096):
        img = Image.open('real_limb.png')
        sun_radius = all_wavelength_dict[self.wavelength]['sun_radius']
        val = int((sun_radius / all_wavelength_dict[self.wavelength][
            'pixel_arcsec_x']) * 2)  # divide 8 to lower resolution and multiply 2 cuz diamter = radius*2
        limb_image_data = img.resize((val, val))
        limb_image_data = np.rot90(limb_image_data)
        pad = int((full_image_val - val) / 2)
        self.limb_image_data = np.pad(limb_image_data, ((pad, pad), (pad, pad), (0, 0)), 'constant',
                                      constant_values=np.nan)

    def draw_cropped_limb(self, bruh=300):
        self.x = int(self.x)
        self.y = int(self.y)
        cropped_limb_image_data = self.limb_image_data[self.x - bruh:self.x + bruh, self.y - bruh:self.y + bruh]
        self.limb_data = cropped_limb_image_data
        self.limb_x = 0
        self.limb_y = 0

    def draw_1600_contours_on_widget(self, contour_display_list):
        self.clear_1600_contours()

        for polygon in contour_display_list:
            self.ax2D.addItem(polygon)
            self.contour_display_filter.append(polygon)

    def clear_1600_contours(self):
        try:
            for linepath in self.contour_display_filter:
                self.ax2D.removeItem(linepath)


        except Exception as e:
            print(e)
            print("exceptionnn")

        self.contour_display_filter = []

    def draw_HMI_polarities_on_widget(self, polarity_display_list):

        self.clear_HMI_contours()
        polygon_list, negative_polygon_list = polarity_display_list

        for polygon in polygon_list:
            self.ax2D.addItem(polygon)
            self.polarity_display_filter.append(polygon)

        for negative_polygon in negative_polygon_list:
            self.ax2D.addItem(negative_polygon)
            self.polarity_display_filter.append(negative_polygon)

    def clear_HMI_contours(self):
        try:
            for linepath in self.polarity_display_filter:
                self.ax2D.removeItem(linepath)

        except Exception as e:
            print(e)
            print("exceptionnn")

        self.polarity_display_filter = []

    def initialize_cropped_image_axis(self):
        self.view.setLimits(xMin=0, xMax=600, yMin=0, yMax=600)
        self.update_axis_dict(peak_x, peak_y)

    def initialize_contour_menu(self):
        self.contour_hide_or_show = False
        contours_button = self.view.menu.addAction('Show 1600 Ribbons')
        contours_button.setCheckable(True)
        contours_button.triggered.connect(self.contour_1600_display_func)

    def initialize_hmi_center_menu(self):
        self.hmi_center_hide_or_show = False
        HMI_center_button = self.view.menu.addAction('Show HMI Center')
        HMI_center_button.setCheckable(True)
        HMI_center_button.triggered.connect(self.show_hmi_center_func)

    def show_hmi_center(self):
        x, y = self.full_HMI_image.get_HMI_centroid_coords()
        x = x + 300 - self.x
        y = y + 300 - self.y
        self.hmi_center_mark1 = LineItem(QtCore.QLineF(x - 20, y - 20, x + 20, y + 20))
        self.hmi_center_mark2 = LineItem(QtCore.QLineF(x - 20, y + 20, x + 20, y - 20))
        pen = QPen()  # creates a default pen
        pen.setWidth(3)
        pen.setBrush(Qt.green)
        self.hmi_center_mark1.setPen(pen)
        self.hmi_center_mark2.setPen(pen)
        self.ax2D.addItem(self.hmi_center_mark1)
        self.ax2D.addItem(self.hmi_center_mark2)
        self.hmi_center_hide_or_show = True

    def hide_hmi_center(self):
        self.ax2D.removeItem(self.hmi_center_mark1)
        self.ax2D.removeItem(self.hmi_center_mark2)
        self.hmi_center_hide_or_show = False

    def show_hmi_center_func(self):
        if self.hmi_center_hide_or_show:
            self.hide_hmi_center()
        else:
            self.show_hmi_center()

    def initialize_polarities_menu(self):
        self.polarities_hide_or_show = False

        polarities_button = self.view.menu.addAction('Show HMI Polarities')
        polarities_button.setCheckable(True)
        polarities_button.triggered.connect(self.polarities_HMI_display_func)

    def contour_1600_display_func(self):
        if self.contour_hide_or_show == True:
            self.contour_hide_or_show = False
            self.clear_1600_contours()

        else:
            self.contour_hide_or_show = True

    def polarities_HMI_display_func(self):
        if self.polarities_hide_or_show == True:
            self.clear_HMI_contours()
            self.polarities_hide_or_show = False
        else:
            self.polarities_hide_or_show = True

    def create_contour_image_data(self, x_click, y_click):
        image_list = all_wavelength_dict[self.wavelength]['full_raw_image_data']
        new_image_list = []
        bruh = 300
        x_click = int(x_click)
        y_click = int(y_click)
        for image in image_list:
            cropped_image_data = image[x_click - bruh:x_click + bruh, y_click - bruh:y_click + bruh]
            new_image_list.append(cropped_image_data)
        all_wavelength_dict[self.wavelength]['raw_cropped_image_data'] = new_image_list
        return new_image_list


class Cropped_RHESSI_Map(Cropped_Maps):
    def __init__(self, win, window, wavelength, size, full_HMI_image):
        super().__init__(win, window, wavelength, size, full_HMI_image)
        self.change_cropped_RHESSI_map_axis(RHESSI_center_x, RHESSI_center_y)

        self.time_label.setText(all_wavelength_dict['RHESSI']['time_string_list'][0])
        self.wavelength_label.setText("RHESSI")
        self.putLabel(self.time_label, self.ax2D, x=3, y=8)
        self.put_wavelength_label(self.wavelength_label, self.ax2D, x=3, y=70)

        self.x = RHESSI_center_x
        self.y = RHESSI_center_y
        self.update_rhessi_axis_dict(RHESSI_center_x, RHESSI_center_y)
        self.draw_RHESSI_limb()
        self.initialize_cropped_limb(full_image_val=512)
        self.draw_cropped_limb(bruh=37)

    def update_rhessi_axis_dict(self, ix, iy):
        self.view.setLimits(xMin=0, xMax=74, yMin=0, yMax=74)
        square_side_length = 37
        small_x = ix - square_side_length
        small_y = iy - square_side_length
        big_x = ix + square_side_length
        big_y = iy + square_side_length

        helioprojective_small_x = RHESSI_pixel_to_helioprojective(small_x, self.sun_center_x, self.pixel_arcsec_x)
        helioprojective_small_y = RHESSI_pixel_to_helioprojective(small_y, self.sun_center_x, self.pixel_arcsec_y)

        rounded_small_x = roundup(helioprojective_small_x)
        rounded_small_y = rounddown(helioprojective_small_y)

        first_x = RHESSI_helioprojective_to_pixel(rounded_small_x, self.sun_center_x, self.pixel_arcsec_x, 0) - small_x
        second_x = RHESSI_helioprojective_to_pixel(rounded_small_x, self.sun_center_x, self.pixel_arcsec_x,
                                                   100) - small_x
        third_x = RHESSI_helioprojective_to_pixel(rounded_small_x, self.sun_center_x, self.pixel_arcsec_x,
                                                  200) - small_x
        first_y = RHESSI_helioprojective_to_pixel(rounded_small_y, self.sun_center_y, self.pixel_arcsec_y, 0) - small_y
        second_y = RHESSI_helioprojective_to_pixel(rounded_small_y, self.sun_center_y, self.pixel_arcsec_y,
                                                   100) - small_y
        third_y = RHESSI_helioprojective_to_pixel(rounded_small_y, self.sun_center_y, self.pixel_arcsec_y,
                                                  200) - small_y

        if (rounded_small_x - helioprojective_small_x) < 59:
            fourth_x = x_helioprojective_to_pixel(rounded_small_x, 300) - small_x
            x_labels = [(first_x, str(rounded_small_x) + '"'), (second_x, str(rounded_small_x + 100) + '"'),
                        (third_x, str(rounded_small_x + 200) + '"'), (fourth_x, str(rounded_small_x + 300) + '"')]
        else:
            x_labels = [(first_x, str(rounded_small_x) + '"'), (second_x, str(rounded_small_x + 100) + '"'),
                        (third_x, str(rounded_small_x + 200) + '"')]

        if (rounded_small_y - helioprojective_small_y) < 1000:
            fourth_y = y_helioprojective_to_pixel(rounded_small_y, 300) - small_y
            y_labels = [(first_y, str(rounded_small_y) + '"'), (second_y, str(rounded_small_y + 100) + '"'),
                        (third_y, str(rounded_small_y + 200) + '"'), (fourth_y, str(rounded_small_y + 300) + '"')]

        else:
            y_labels = [(first_y, str(rounded_small_y) + '"'), (second_y, str(rounded_small_y + 100) + '"'),
                        (third_y, str(rounded_small_y + 200) + '"')]

        self.axis_left.setTicks([y_labels])
        self.axis_bottom.setTicks([x_labels])


class Cropped_171_Map(Cropped_Maps):
    def __init__(self, win, window, wavelength, size, full_HMI_image):
        super().__init__(win, window, wavelength, size, full_HMI_image)
        self.initialize_show_center_button()


class Cropped_1600_Map(Cropped_Maps):
    def __init__(self, win, window, wavelength, size, full_HMI_image):
        super().__init__(win, window, wavelength, size, full_HMI_image)
        self.initialize_1600_sunpy_map()
        self.create_contours_list(peak_x, peak_y)



    def initialize_1600_sunpy_map(self):
        self.cea_1600_list = []
        for i in range(0,int(len(all_wavelength_dict[1600]['full_raw_image_data']))):

            data = all_wavelength_dict[1600]['full_raw_image_data'][i]
            header = all_wavelength_dict[1600]['header_list'][i]
            data = np.rot90(data)
            data = np.flipud(data)
            map_1600 = sunpy.map.Map(data, header)

            x, y = self.full_HMI_image.get_hmi_center()
            #map_1600_cea, _ = HMI_to_CEA_func(map_1600, 512, 1024)

            #HMI_positive_array = copy.deepcopy(self.full_HMI_image.hmi_array)
            """
            top_right = SkyCoord((x + 10) * u.deg, (y + 5) * u.deg,
                                 frame=map_1600_cea.coordinate_frame)
            bottom_left = SkyCoord((x - 10) * u.deg,
                                   (y - 5) * u.deg, frame=map_1600_cea.coordinate_frame)

            submap_1600 = (map_1600_cea.submap(bottom_left, top_right=top_right))
            """

            self.cea_1600_list.append(map_1600)


    def calculate_flux_weighted_centroid(self):
        self.ribbon_weighted_centroids_dict = {}
        for i in range(0,len(self.image_list)):
            start = lmao.time()
            array = copy.deepcopy(self.image_list[i])
            x,y = HMI_weighted_centering_func(array,'greater',self.max_clip_val)
            x = 600-x
            y = 600-y
            print(type(x))
            if str(x) == "nan" or str(y) == "nan":
                x,y,helioprojective_x,helioprojective_y,heliographic_coords = "N/A","N/A","N/A","N/A","N/A"
            else:
                x = self.x - 300 + x
                y = self.y - 300 + y
                helioprojective_x = x_pixel_to_helioprojective(x, self.sun_center_x, self.pixel_arcsec_x)
                helioprojective_y = y_pixel_to_helioprojective(y, self.sun_center_x, self.pixel_arcsec_y)

                if helioprojective_x > 1100 or helioprojective_x < -1100 or helioprojective_y > 1100 or helioprojective_y < -1100:
                    helioprojective_x = 0
                    helioprojective_y = 0
                print("index")
                print(i)
                print(self.max_clip_val)
                print(x, y)
                heliographic_coords = helioprojective_to_heliographic(helioprojective_x, helioprojective_y,
                                                                      self.cea_1600_list[i])
                end = lmao.time()
            self.ribbon_weighted_centroids_dict[i] = {}
            self.ribbon_weighted_centroids_dict[i]['pixel'] = x, y
            self.ribbon_weighted_centroids_dict[i]['helioprojective'] = helioprojective_x, helioprojective_y
            self.ribbon_weighted_centroids_dict[i]['heliographic'] = heliographic_coords

    def create_contours_list(self, x, y):
        x, y = int(x), int(y)
        self.image_list = self.create_contour_image_data(x, y)
        self.calculate_flux_weighted_centroid()
        self.aia_1600_contour_list = []
        for data in self.image_list:
            data = np.rot90(data)
            data = np.flipud(data)
            fig1, ax2 = plt.subplots()

            b = ax2.contour(data, [int(self.max_clip_val/2), 1000000], origin='lower', colors=['red'], linestyles='solid')

            p1 = b.collections[0].get_paths()
            p1.sort(key=len, reverse=True)

            contour_list = []

            try:
                for i in range(0, 5):
                    negative_vertices = p1[i].vertices
                    contour_list.append(negative_vertices)

            except IndexError:
                print("less than 5 contours")
                pass

            contour_list = convert_to_Qpoint_list(contour_list, 1)

            self.aia_1600_contour_list.append(contour_list)

            plt.close()

        self.convert_to_polygon_1600(self.aia_1600_contour_list)
        self.image_list = asinh_stretch_images(ct.aia_color_table(1600 * u.angstrom), 0, 1000, self.image_list)

    def convert_to_polygon_1600(self, contour_list):
        self.contour_display_list = []
        self.contour_display_list1 = []
        self.contour_display_list2 = []
        self.contour_display_list3 = []
        lists = [self.contour_display_list, self.contour_display_list1, self.contour_display_list2,
                 self.contour_display_list3]
        for list1 in lists:
            for contours in contour_list:
                polygon_list = []
                for linepath in contours:
                    polygon = create_polygon(linepath, QtCore.Qt.green)
                    polygon_list.append(polygon)

                list1.append(polygon_list)


class Cropped_Helioprojective_HMI_Map(Cropped_Maps):
    def __init__(self, win, window, wavelength, size, full_HMI_image):
        super().__init__(win, window, wavelength, size, full_HMI_image)
        self.create_contours_list(peak_x, peak_y)

    def create_contours_list(self, x, y):
        self.x = x
        self.y = y
        self.image_list = self.clip_image_list(self.create_contour_image_data(x, y), -50, 50)
        self.HMI_polarity_list = []

        for data in self.image_list:
            data = np.rot90(data)
            data = np.flipud(data)
            fig1, ax2 = plt.subplots()
            b = ax2.contour(data, [-50, -49], origin='lower', colors=['red'], linestyles='solid')
            a = ax2.contour(data, [49, 50], origin='lower', colors=['blue'], linestyles='solid')

            p1 = a.collections[0].get_paths()
            p1.sort(key=len, reverse=True)
            p2 = b.collections[0].get_paths()
            p2.sort(key=len, reverse=True)

            positive_contour_list = []
            negative_contour_list = []

            try:
                for i in range(0, 10):
                    negative_vertices = p1[i].vertices
                    negative_contour_list.append(negative_vertices)
                    positive_vertices = p2[i].vertices
                    positive_contour_list.append(positive_vertices)
            except IndexError:
                print("less than 10 contours")
                pass

            cropped_positive_list = convert_to_Qpoint_list(positive_contour_list, 1)
            cropped_negative_list = convert_to_Qpoint_list(negative_contour_list, 1)

            self.HMI_polarity_list.append([cropped_positive_list, cropped_negative_list])

            plt.close()

        self.convert_to_polygon_HMI(self.HMI_polarity_list)

    def convert_to_polygon_HMI(self, contour_list):
        self.polarity_display_list = []
        self.polarity_display_list1 = []
        self.polarity_display_list2 = []
        self.polarity_display_list3 = []

        lists = [self.polarity_display_list, self.polarity_display_list1, self.polarity_display_list2,
                 self.polarity_display_list3]
        for list1 in lists:
            for contours in contour_list:
                positive_list = []
                negative_list = []

                cropped_positive_list, cropped_negative_list = contours

                for linepath in cropped_positive_list:
                    polygon = create_polygon(linepath, QtCore.Qt.blue)
                    positive_list.append(polygon)

                for linepath in cropped_negative_list:
                    polygon = create_polygon(linepath, QtCore.Qt.red)
                    negative_list.append(polygon)

                list1.append([positive_list, negative_list])


class Full_HMI_Image():
    def __init__(self, ui):
        self.ui = ui
        self.initialize_hmi_text_boxes()
        self.aia_map = aia_map_object
        self.HMI_full_image_func()

    def initialize_hmi_text_boxes(self):
        self.user_hmi_x_coord = self.ui.user_hmi_x_coord
        self.user_hmi_y_coord = self.ui.user_hmi_y_coord

    def get_HMI_centroid_coords(self):
        heliographic_coord = self.swap_submap.pixel_to_world(self.x_center_hmi * u.pix, self.y_center_hmi * u.pix)

        helio_coord = heliographic_coord.transform_to(self.aia_map.coordinate_frame)

        pixel_coord = self.aia_map.world_to_pixel(helio_coord)

        x_coord = int(float(str(pixel_coord[0])[:-3]))
        y_coord = int(float(str(pixel_coord[1])[:-3]))
        return x_coord, y_coord

    def get_hmi_center(self):
        return self.x_hmi_coord, self.y_hmi_coord

    def HMI_full_image_func(self):

        scene = self.ui.scene

        self.m_hmi_cea, x_negative1, x_positive1, y_negative1, y_positive1, hmi_center_x_pixel, hmi_center_y_pixel, full_data, self.m_hmi_cea_copy, self.wcs_axes, x_center_heliographical_coord, y_center_heliographical_coord = HMI_object
        fig = Figure(figsize=(6.25, 3.515625))

        fig.subplots_adjust(bottom=0.2)
        self.ax11 = fig.gca(projection=self.m_hmi_cea)

        self.ax11.set_xlabel('Heliographic Latitude', color='white')
        self.ax11.set_ylabel('Longitude', color='white')

        fig.patch.set_facecolor('black')
        data = np.clip(full_data, -50, 50)
        self.ax11.tick_params(axis='x', colors='white')
        self.ax11.tick_params(axis='y', colors='white')
        c = self.ax11.imshow(data, cmap='gray')

        self.canvas = FigureCanvas(fig)
        proxy_widget = scene.addWidget(self.canvas)


        self.user_hmi_x_coord.setText(str(round(x_center_heliographical_coord, 4)))
        self.user_hmi_y_coord.setText(str(round(y_center_heliographical_coord, 4)))
        self.x_hmi_coord = float(x_center_heliographical_coord)
        self.y_hmi_coord = float(y_center_heliographical_coord)
        self.HMI_cropped_image_func(x_center_heliographical_coord, y_center_heliographical_coord)

        fig.canvas.mpl_connect('button_press_event', self.HMI_onclick)

    def HMI_onclick(self, event):
        global user_clicked_map_bool
        if (event.inaxes is not None):
            ix, iy = int(event.xdata), int(event.ydata)

            heliographic_coord = self.m_hmi_cea_copy.pixel_to_world(ix * u.pix, iy * u.pix)

            EW_lon, NS_lat = heliographic_coord.to_string('decimal').split(' ')

            EW_lon = float(EW_lon)
            NS_lat = float(NS_lat)

            self.x_hmi_coord = EW_lon
            self.y_hmi_coord = NS_lat
            self.user_hmi_x_coord.setText(str(self.x_hmi_coord))
            self.user_hmi_y_coord.setText(str(self.y_hmi_coord))
            self.HMI_cropped_image_func(EW_lon, NS_lat)

    def calculate_area(self, cdelt1, cdelt2, num_pixels):
        dx = cdelt1 * (0.5 / 180) * 3600 * 727
        dy = cdelt2 * (0.5 / 180) * 3600 * 727
        da = dx * dy  # x times y is area of pixel
        area_ar = num_pixels * da  # area in km^2
        area_calc = area_ar / 1000000

        return area_calc

    def calculate_HMI_area(self, submap):
        lv1 = 49  # Standard deviation of the B distribution on the SHARP region
        submap_data = np.clip(submap.data, -50, 50)
        num_pixels = np.where(abs(submap_data) >= lv1)[0].size

        meta = submap.meta
        cdelt1 = float(meta['cdelt1'])
        cdelt2 = float(meta['cdelt2'])
        area_calc = self.calculate_area(cdelt1, cdelt2, num_pixels)
        return area_calc, num_pixels

    def calculate_ribbon_area(self, HMI_positive_array, submap, max_clip_val):
        max_clip_val = int(max_clip_val/2)
        submap_1600_data = submap.data
        self.ax12.contour(submap_1600_data, [max_clip_val, 999999999999], origin='lower', colors=['green'],
                          linestyles='solid')
        self.canvas2.draw_idle()
        HMI_negative_array = copy.deepcopy(HMI_positive_array)
        HMI_positive_array[HMI_positive_array < 49] = 0
        HMI_positive_array[HMI_positive_array >= 49] = 1
        HMI_negative_array[HMI_negative_array > -49] = 0
        HMI_negative_array[HMI_negative_array <= -49] = 1

        submap_1600_data[submap_1600_data < max_clip_val] = -1
        submap_1600_data[submap_1600_data >= max_clip_val] = 1

        positive_counter = 0

        minimum_x = min([HMI_negative_array.shape[0], HMI_positive_array.shape[0], submap_1600_data.shape[0]])
        minimum_y = max([HMI_negative_array.shape[0], HMI_positive_array.shape[0], submap_1600_data.shape[0]])

        for i in range(0, minimum_x):
            for j in range(0, minimum_y):
                if HMI_positive_array[i][j] == submap_1600_data[i][j]:
                    positive_counter += 1

        negative_counter = 0
        for i in range(0, minimum_x):
            for j in range(0, minimum_y):
                if HMI_negative_array[i][j] == submap_1600_data[i][j]:
                    negative_counter += 1
        meta = submap.meta
        cdelt1 = float(meta['cdelt1'])
        cdelt2 = float(meta['cdelt2'])
        positive_ribbons_area = self.calculate_area(cdelt1, cdelt2, positive_counter)
        negative_ribbons_area = self.calculate_area(cdelt1, cdelt2, negative_counter)

        return positive_ribbons_area, negative_ribbons_area

    def change_hmi_view_func(self):
        user_hmi_x = float(self.ui.user_hmi_x_coord.toPlainText())
        user_hmi_y = float(self.ui.user_hmi_y_coord.toPlainText())

        x, y = self.HMI_cropped_image_func(user_hmi_x, user_hmi_y)

        x, y = self.convert_hmi_centroid_to_aia_pixel(x, y)
        return x, y

    def convert_hmi_centroid_to_aia_pixel(self, x, y):
        heliographic_coord = self.swap_submap.pixel_to_world(x * u.pix, y * u.pix)
        helio_coord = heliographic_coord.transform_to(self.aia_map.coordinate_frame)
        pixel_coord = self.aia_map.world_to_pixel(helio_coord)
        x_coord = int(float(str(pixel_coord[0])[:-3]))
        y_coord = int(float(str(pixel_coord[1])[:-3]))
        return x_coord, y_coord

    def HMI_cropped_image_func(self, x_center_heliographical_coord, y_center_heliographical_coord):
        try:
            self.ax11.patches.pop()
        except IndexError:
            pass

        top_right = SkyCoord((x_center_heliographical_coord + 10) * u.deg, (y_center_heliographical_coord + 5) * u.deg,
                             frame=self.m_hmi_cea.coordinate_frame)
        bottom_left = SkyCoord((x_center_heliographical_coord - 10) * u.deg,
                               (y_center_heliographical_coord - 5) * u.deg, frame=self.m_hmi_cea.coordinate_frame)

        self.swap_submap = self.m_hmi_cea_copy.submap(bottom_left, top_right=top_right)
        pd.to_pickle(self.swap_submap,"swap_submap.pickle")
        area_number, num_pixels = self.calculate_HMI_area(self.swap_submap)
        flux_sum = abs(self.swap_submap.data).sum()
        HMI_flux = int(area_number) * flux_sum/ num_pixels * 10e10*1000000

        self.ui.HMI_area_value.setText(str(int(area_number)) + "Mega-meters²")
        self.ui.HMI_flux_value.setText(str("{:.2e}".format(int(HMI_flux))) + "Mx")

        scene2 = self.ui.scene2

        try:
            self.ax12.clear()
            self.fig2.clf()

        except AttributeError:
            pass
        self.fig2 = Figure(figsize=(6.25, 3.515625))
        self.fig2.subplots_adjust(bottom=0.2)
        self.fig2.patch.set_facecolor('black')

        self.canvas2 = FigureCanvas(self.fig2)
        scene2.addWidget(self.canvas2)
        self.ax12 = self.fig2.gca(projection=self.swap_submap)

        data = np.clip(self.swap_submap.data, -50, 50)
        self.ax12.imshow(data, cmap='gray', aspect='auto')

        positive_array = self.swap_submap.data
        negative_array = copy.deepcopy(positive_array)
        self.hmi_array = negative_array

        x_positive, y_positive = HMI_weighted_centering_func(positive_array, 'greater', 50)

        x_negative, y_negative = HMI_weighted_centering_func(negative_array, 'less', 50)
        self.x_center_hmi = (x_positive + x_negative) / 2
        self.y_center_hmi = (y_positive + y_negative) / 2

        self.ax12.plot(x_positive, y_positive, 'wx', color='red', markersize=20)
        self.ax12.plot(x_negative, y_negative, 'wx', color='blue', markersize=20)
        self.ax12.plot(self.x_center_hmi, self.y_center_hmi, 'wx', color='yellow', markersize=20)

        self.ax12.set_xlabel('Heliographic Latitude', color='white')
        self.ax12.set_ylabel('Longtitude', color='white')
        self.ax12.tick_params(axis='x', colors='white')
        self.ax12.tick_params(axis='y', colors='white')



        b = self.m_hmi_cea_copy.world_to_pixel(bottom_left)
        x_coord = float(str(b[0])[:-3])
        y_coord = float(str(b[1])[:-3])
        rect = patches.Rectangle((x_coord - 50, y_coord),
                                 300,
                                 150, linewidth=1.5, edgecolor='g', facecolor='none')

        self.ax11.add_patch(rect)
        self.ax12.contour(data, [-50, -49], origin='lower', colors=['blue'], linestyles='solid')
        self.ax12.contour(data, [49, 50], origin='lower', colors=['red'], linestyles='solid')

        self.canvas.draw_idle()
        self.canvas2.draw_idle()

        top_right = SkyCoord((x_center_heliographical_coord + 10) * u.deg, (y_center_heliographical_coord + 5) * u.deg,
                             frame=self.m_hmi_cea_copy.coordinate_frame)
        bottom_left = SkyCoord((x_center_heliographical_coord - 10) * u.deg,
                               (y_center_heliographical_coord - 5) * u.deg, frame=self.m_hmi_cea_copy.coordinate_frame)

        self.swap_submap = self.m_hmi_cea_copy.submap(bottom_left, top_right=top_right)

        return self.x_center_hmi, self.y_center_hmi


class Graph():
    def __init__(self, date_list, intensity_list, vb_1, vb_2, clickbox, color):
        self.intensity_list = intensity_list
        self.vb_1 = vb_1
        self.vb_2 = vb_2
        self.clickbox = clickbox
        self.color = color
        self.date_list = []
        for date1 in date_list:
            time = (date1 - datetime(1969, 12, 31, 17, 0, 0)).total_seconds()
            self.date_list.append(time)

        self.create_graph()


        self.graph.setLogMode(False, False)  # y is logarithmic
        self.clickbox.setChecked(True)
        self.clickbox.stateChanged.connect(self.clickBox_func)
        self.vb_2.addItem(self.graph)
        self.vb_2.setXLink(vb_1)

    def create_graph(self):
        self.graph = pg.PlotDataItem(self.date_list, self.intensity_list,
                                     pen=pg.mkPen(self.color, width=1))

    def update_data(self):
        self.graph.setData(self.date_list, self.intensity_list)

    def clickBox_func(self, state):
        if state == QtCore.Qt.Checked:
            self.graph.setVisible(True)
        else:
            self.graph.setVisible(False)

    def uncheck_func(self):
        self.graph.setVisible(False)
        self.clickbox.setChecked(False)

    def check_func(self):
        self.graph.setVisible(True)
        self.clickbox.setChecked(True)


class MainWindow():
    def __init__(self):
        self.painter = QPainter()

        self.i = 0
        self.frame_timer = 20

        self.full_HMI_map_time_label = pg.TextItem('helo', **{'color': '#ffffff'})
        self.cropped_HMI_map_time_label = pg.TextItem('helo', **{'color': '#ffffff'})
        self.full_RHESSI_time_label = pg.TextItem('helo', **{'color': '#ffffff'})

        self.main_window = QMainWindow()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self.main_window)
        self.create_list_slice_dates()
        self.initialize_radio_buttons()
        self.initialize_metadata()
        self.initialize_time_list()
        self.goes_euv_plot()
        self.initialize_widgets()
        self.coordinate_system = 'Helioprojective'
        self.full_1600_image_widget_func()
        self.movie_controls()
        self.update_RHESSI_centering_label()
        self.initialize_aia_centering_label(peak_x, peak_y)
        self.create_list_slice_dates()
        self.run_widgets_func()

    def pickle_func(self):
        measurements_data = pd.read_csv('measurements.csv',index_col='index')
        euv_dict = {}
        for wavelength1 in wavelength_list:
            dict1 = {}
            intensity_list = self.EUV_graph_dict[wavelength1].intensity_list
            time_list = all_wavelength_dict[wavelength1]['time_string_list']
            for i in range(0,len(intensity_list)):
                dict1[time_list[i]] = intensity_list[i]

            euv_dict[wavelength1] = dict1
        measurement_dict = {'GOES Start Time':self.ui.start_time_value.toPlainText(), 'GOES Peak Time':self.ui.peak_time_value.toPlainText(), 'GOES End Time':self.ui.end_time_value.toPlainText(), 'Flare Class': self.ui.flare_class_value.toPlainText(), 'CME': self.ui.CME_boolean.toPlainText(), 'CME Speed': self.ui.CME_speed_value.toPlainText(),
                            'CME Angle': self.ui.CME_angle_value.toPlainText(), "+ Ribbon Area": self.ui.positive_ribbon_area_value.toPlainText(),"- Ribbon Area": self.ui.negative_ribbon_area_value.toPlainText(),'HMI Flux': self.ui.HMI_flux_value.toPlainText(), 'Active Region #': self.ui.active_region_number.toPlainText(), 'Hale Class': self.ui.hale_class_value.toPlainText(),
                            'EUV Data': euv_dict, 'Ribbon Centroid per frame': self.cropped_1600_image_widget.ribbon_weighted_centroids_dict, 'HMI Ribbon Centroid per frame': self.HMI_ribbons_centroid_dict}
        measurements_data.loc[len(measurements_data)] = measurement_dict
        measurements_data.to_csv('measurements.csv')


    def calculate_HMI_ribbons_centroid(self):
        self.HMI_ribbons_centroid_dict = {}
        cropped_image_list_HMI = copy.deepcopy(self.cropped_helioprojective_HMI_image_widget.image_list)
        cropped_image_list_1600 = copy.deepcopy(all_wavelength_dict[1600]['raw_cropped_image_data'])

        for i in range(0,len(cropped_image_list_HMI)):

            HMI_positive_array = cropped_image_list_HMI[i]
            HMI_negative_array = copy.deepcopy(HMI_positive_array)

            submap_1600_data = cropped_image_list_1600[i]


            HMI_positive_array[HMI_positive_array < 49] = 0
            HMI_positive_array[HMI_positive_array >= 49] = 1
            HMI_negative_array[HMI_negative_array > -49] = 0
            HMI_negative_array[HMI_negative_array <= -49] = 1



            submap_1600_data = np.where(submap_1600_data < all_wavelength_dict[1600]['max_clip_val'], -1, submap_1600_data)
            submap_1600_data = np.where(submap_1600_data >= all_wavelength_dict[1600]['max_clip_val'], 1, submap_1600_data)
            print("all wavelength")
            print(all_wavelength_dict[1600]['max_clip_val'])

            minimum_x = min([HMI_negative_array.shape[0], HMI_positive_array.shape[0], submap_1600_data.shape[0]])
            minimum_y = max([HMI_negative_array.shape[0], HMI_positive_array.shape[0], submap_1600_data.shape[0]])

            positive_overlapping_array = np.zeros([600,600])
            negative_overlapping_array = np.zeros([600, 600])

            p = 0
            for k in range(0, minimum_x):
                for j in range(0, minimum_y):
                    if HMI_positive_array[k][j] == submap_1600_data[k][j]:
                        positive_overlapping_array[k][j] = 1
                        p+= 1
            print(p)
            print("bruh")
            for k in range(0, minimum_x):
                for j in range(0, minimum_y):
                    if HMI_negative_array[k][j] == submap_1600_data[k][j]:
                        negative_overlapping_array[k][j] = 1


            positive_y, positive_x = ndimage.measurements.center_of_mass(positive_overlapping_array)
            negative_y, negative_x = ndimage.measurements.center_of_mass(negative_overlapping_array)
            if str(positive_y) == "nan" or str(negative_x) == "nan":
                x, y, helioprojective_x, helioprojective_y, heliographic_coords = "N/A", "N/A", "N/A", "N/A", "N/A"
            else:

                x, y = (positive_x + negative_x) / 2, (negative_x + negative_y) / 2
                x = self.cropped_1600_image_widget.x - 300 + (600 - x)
                y = self.cropped_1600_image_widget.y - 300 + (600 - y)
                helioprojective_x = x_pixel_to_helioprojective(x, self.cropped_1600_image_widget.sun_center_x,
                                                               self.cropped_1600_image_widget.pixel_arcsec_x)
                helioprojective_y = y_pixel_to_helioprojective(y, self.cropped_1600_image_widget.sun_center_x,
                                                               self.cropped_1600_image_widget.pixel_arcsec_y)
                print("index")
                print(i)

                if helioprojective_x > 1100 or helioprojective_x < -1100 or helioprojective_y > 1100 or helioprojective_y < -1100:
                    helioprojective_x = 0
                    helioprojective_y = 0

                heliographic_coords = helioprojective_to_heliographic(helioprojective_x, helioprojective_y,
                                                                      self.cropped_1600_image_widget.cea_1600_list[i])

            self.HMI_ribbons_centroid_dict[i] = {}
            self.HMI_ribbons_centroid_dict[i]['pixel'] = x, y
            self.HMI_ribbons_centroid_dict[i]['helioprojective'] = helioprojective_x, helioprojective_y
            self.HMI_ribbons_centroid_dict[i]['heliographic'] = heliographic_coords

    def initialize_metadata(self):
        self.ui.start_time_value.setText(flare_metadata['GOES Start Time'][:19])
        self.ui.peak_time_value.setText(flare_metadata['GOES Peak Time'][:19])
        self.ui.end_time_value.setText(flare_metadata['GOES End Time'][:19])
        self.ui.CME_angle_value.setText(str(flare_metadata['Angle']))
        self.ui.CME_speed_value.setText(str(flare_metadata['Speed']))
        self.ui.CME_boolean.setText(str(flare_metadata['CME']))
        self.ui.flare_class_value.setText((flare_metadata['Flare Class']))

    def initialize_widgets(self):
        self.full_rhessi_widget = Full_RHESSI_Widget(win=self.ui.RHESSI_image_widget, window=self.main_window,
                                                     size='full', wavelength='RHESSI')

        self.full_HMI_image = Full_HMI_Image(self.ui)
        self.cropped_aia_map_widget = Cropped_171_Map(win=self.ui.cropped_aia_map_widget, window=self.main_window,
                                                      size='cropped', wavelength=171,
                                                      full_HMI_image=self.full_HMI_image)

        self.cropped_rhessi_map_widget = Cropped_RHESSI_Map(win=self.ui.RHESSI_cropped_image_widget,
                                                            window=self.main_window,
                                                            size='cropped', wavelength='RHESSI',
                                                            full_HMI_image=self.full_HMI_image)

        self.cropped_1600_image_widget = Cropped_1600_Map(win=self.ui.cropped_1600_image_widget,
                                                          window=self.main_window, size='cropped', wavelength=1600,
                                                          full_HMI_image=self.full_HMI_image)

        self.cropped_helioprojective_HMI_image_widget = Cropped_Helioprojective_HMI_Map(
            win=self.ui.cropped_helioprojective_HMI_image_widget,
            window=self.main_window, size='cropped', wavelength='HMI', full_HMI_image=self.full_HMI_image)

        self.full_1600_image_widget = Other_AIA_Widgets(win=self.ui.full_1600_image_widget, window=self.main_window,
                                                        size='full', wavelength=1600)

        self.HMI_image_widget = Full_Helioprojective_HMI_Widget(win=self.ui.HMI_image_widget, window=self.main_window,
                                                                size='full', wavelength='HMI')

        self.full_aia_map_widget = Full_AIA_Map_Widget(
            wavelength=171, size1='full', window=self.main_window,
            ui_hale_class_value=self.ui.hale_class_value,
            ui_active_region_number=self.ui.active_region_number,
            cropped_map_object=self.cropped_aia_map_widget, full_hmi_map_object=self.full_HMI_image,
            cropped_1600_map_object=self.cropped_1600_image_widget,
            cropped_hmi_map_object=self.cropped_helioprojective_HMI_image_widget,
            full_helioprojective_hmi_map_object=self.HMI_image_widget, full_1600_map_object=self.full_1600_image_widget,graph_list=self.EUV_graph_dict,normalized_graph_list=self.normalized_EUV_graph_dict)

        self.full_aia_map_widget.setParent(self.ui.tab_1)
        self.full_aia_map_widget.setGeometry(QtCore.QRect(0, 383, 347, 347))
        self.full_aia_map_widget.setObjectName("full_aia_map_widget")

        self.aia_widget_0 = Other_AIA_Widgets(win=self.ui.aia_widget_0, window=self.main_window, size='full',
                                              wavelength=94)
        self.aia_widget_3 = Other_AIA_Widgets(win=self.ui.aia_widget_3, window=self.main_window, size='full',
                                              wavelength=211)
        self.aia_widget_4 = Other_AIA_Widgets(win=self.ui.aia_widget_4, window=self.main_window, size='full',
                                              wavelength=304)


        self.initialize_wavelength_and_state_values()

        self.widget_list = [self.cropped_aia_map_widget, self.full_aia_map_widget, self.aia_widget_0, self.aia_widget_3, self.aia_widget_4,
                            self.cropped_1600_image_widget, self.full_1600_image_widget, self.HMI_image_widget,
                            self.cropped_helioprojective_HMI_image_widget, self.cropped_rhessi_map_widget]

        self.full_widget_list = [self.full_aia_map_widget, self.aia_widget_0, self.aia_widget_3, self.aia_widget_4,
                                 self.full_1600_image_widget, self.HMI_image_widget]

        self.cropped_widget_list = [self.cropped_aia_map_widget, self.cropped_1600_image_widget,
                                    self.cropped_helioprojective_HMI_image_widget, self.cropped_rhessi_map_widget]

        self.HMI_ribbons_centroid_dict = {}
    def initialize_wavelength_and_state_values(self):
        self.length_of_list = len(self.full_aia_map_widget.time_string_list) - 1

    def get_night_saa_time_range(self, start_time, end_time, csv_file_name):
        night_dataframe = pd.read_csv(csv_file_name)

        night_start_time_list = list(night_dataframe['Start Time'])
        night_end_time_list = list(night_dataframe['End Time'])

        night_start_value = bisect.bisect(night_end_time_list, start_time)

        night_end_value = bisect.bisect(night_start_time_list, end_time) - 1

        import math
        night_index_list = list(range(math.ceil(night_start_value), math.floor(night_end_value) + 1))

        final_index_list = []

        for number in night_index_list:
            start, end = night_dataframe.iloc[number].to_list()
            start = to_datetime_object(start, "%Y-%m-%d %H:%M:%S")
            end = to_datetime_object(end, "%Y-%m-%d %H:%M:%S")
            start = (start - datetime(1969, 12, 31, 17, 0, 0)).total_seconds()
            end = (end - datetime(1969, 12, 31, 17, 0, 0)).total_seconds()
            final_index_list.append([start, end])

        return final_index_list

    def create_list_slice_dates(self):
        self.beginning_time = (
            str(special_beginning_time_str).replace("-", "").replace("T", "").replace(":", "").replace(".", "").replace(
                " ",
                "")[
            :12])
        self.end_time = (
            str(special_end_time_str).replace("-", "").replace("T", "").replace(":", "").replace(".", "").replace(" ",
                                                                                                                  "")[
            :12])

    def initialize_radio_buttons(self):
        print("to pickle")
        self.to_pickle_button = self.ui.to_pickle_button
        self.to_pickle_button.clicked.connect(self.pickle_func)
        self.lograthmic_radio_button = self.ui.log_button_3
        self.lograthmic_radio_button.toggled.connect(self.GOES_scaling_func)
        self.linear_radio_button = self.ui.linear_button_3
        self.linear_radio_button.toggled.connect(self.GOES_scaling_func)
        self.regular_euv_radio_button = self.ui.regular_euv_type_button_3
        self.regular_euv_radio_button.toggled.connect(self.EUV_regular_func)
        self.normalized_euv_radio_button = self.ui.normalized_EUV_type_button_3
        self.normalized_euv_radio_button.toggled.connect(self.EUV_normalized_func)
        self.pixel_button = self.ui.pixel_button_3
        self.pixel_button.toggled.connect(self.coordinate_type_func)
        self.helioprojective_button = self.ui.helioprojective_button_3
        self.helioprojective_button.toggled.connect(self.coordinate_type_func)
        self.heliographical_button = self.ui.heliographical_button_3
        self.heliographical_button.toggled.connect(self.coordinate_type_func)

    def GOES_scaling_func(self, value):
        rbtn = self.main_window.sender()

    def EUV_regular_func(self):
        for graph in self.normalized_EUV_graph_dict.values():
            graph.uncheck_func()
        for graph in self.EUV_graph_dict.values():
            print(graph)
            graph.check_func()
        self.vb_2.disableAutoRange()
        self.vb_2.setLimits(yMin=self.min_intensity_val, yMax=self.max_intensity_val)
        self.vb_2.setGeometry(self.vb_1.sceneBoundingRect())

    def EUV_normalized_func(self):
        for graph in self.normalized_EUV_graph_dict.values():
            graph.check_func()
        for graph in self.EUV_graph_dict.values():
            graph.uncheck_func()
        self.vb_2.disableAutoRange()
        self.vb_2.setLimits(yMin=0, yMax=1)
        self.vb_2.setGeometry(self.vb_1.sceneBoundingRect())

    def coordinate_type_func(self):
        a = self.main_window.sender()
        self.coordinate_system = a.text()
        self.update_box_centering_label()
        self.update_RHESSI_centering_label()

    def update_box_centering_label(self):
        if self.coordinate_system == "Pixel":
            x_coord, y_coord = brightest_pixel_dict['box_position_pixel']
            self.box_position.setText("({},{})".format(str(int(x_coord)), str(int(y_coord))))
        elif self.coordinate_system == "Helioprojective":
            x_coord, y_coord = brightest_pixel_dict['box_position_helioprojective']
            self.box_position.setText("({}'',{}'')".format(str(int(x_coord)), str(int(y_coord))))
        elif self.coordinate_system == "Heliographical":
            coords = brightest_pixel_dict['box_position_heliographical']
            self.box_position.setText(simplify_heliographic_coords(str(coords)))

    def update_RHESSI_centering_label(self):
        self.RHESSI_position_label = self.ui.RHESSI_center_position
        if self.coordinate_system == "Pixel":
            x_coord, y_coord = brightest_pixel_dict['RHESSI_center_pixel']
            self.RHESSI_position_label.setText("({},{})".format(str(int(x_coord)), str(int(y_coord))))
        elif self.coordinate_system == "Helioprojective":
            x_coord, y_coord = brightest_pixel_dict['RHESSI_center_helioprojective']
            self.RHESSI_position_label.setText("({}'',{}'')".format(str(int(x_coord)), str(int(y_coord))))
        elif self.coordinate_system == "Heliographical":
            coords = brightest_pixel_dict['RHESSI_center_heliographical']
            self.RHESSI_position_label.setText(simplify_heliographic_coords(str(coords)))

    def initialize_aia_centering_label(self, x, y):
        self.box_position = self.ui.box_position
        brightest_pixel_dict['box_position_pixel'] = [x, y]
        helioprojective_x = x_pixel_to_helioprojective(x, sun_center_x_list_placeholder[0], pixel_arcsec_x_list_placeholder[0])
        helioprojective_y = y_pixel_to_helioprojective(y, sun_center_y_list_placeholder[0], pixel_arcsec_y_list_placeholder[0])
        brightest_pixel_dict['box_position_helioprojective'] = helioprojective_x, helioprojective_y
        heliographic_coords = helioprojective_to_heliographic(helioprojective_x, helioprojective_y,
                                                              aia_map_object)
        brightest_pixel_dict['box_position_heliographical'] = heliographic_coords.split(',')
        self.update_box_centering_label()
        #self.box_position.setText("({}'',{}'')".format(str(int(helioprojective_x)), str(int(helioprojective_y))))

    def update_centering_labels(self,dict1, label_value, coord_type,frame_number):
        if coord_type == "Pixel":
            x_coord, y_coord = dict1[frame_number]['pixel']
            if not x_coord == "N/A":
                x_coord = str(int(x_coord))
                y_coord = str(int(y_coord))
            label_value.setText("({},{})".format(x_coord,y_coord))
        elif coord_type == "Helioprojective":
            x_coord, y_coord = dict1[frame_number]['helioprojective']
            if not x_coord == "N/A":
                x_coord = str(int(x_coord))
                y_coord = str(int(y_coord))
            label_value.setText("({}'',{}'')".format(x_coord,y_coord))
        elif coord_type == "Heliographical":
            coords = dict1[frame_number]['heliographic']
            label_value.setText(simplify_heliographic_coords(str(coords)))

    def real_time_update_centering_labels(self, frame_number):
        coord_type = self.coordinate_system
        self.update_centering_labels(brightest_pixel_dict,self.ui.center_of_flare,coord_type,frame_number)
        self.update_centering_labels(self.cropped_1600_image_widget.ribbon_weighted_centroids_dict, self.ui.ribbons_1600_position, coord_type, frame_number)
        if not self.HMI_ribbons_centroid_dict == {}:
            self.update_centering_labels(self.HMI_ribbons_centroid_dict, self.ui.HMI_ribbons_position, coord_type, frame_number)

    def putRHESSILabel(self, label, axes):
        label.setPos(QtCore.QPointF(12, 115))
        label.setFont(QFont("Times New Roman", 20))
        axes.addItem(label)

    def GOES_XRSB_clickBox(self, state):
        if state == QtCore.Qt.Checked:
            self.GOES_XRSB_graph.setVisible(True)
        else:
            self.GOES_XRSB_graph.setVisible(False)

    def RHESSI_12_clickBox(self, state):
        if state == QtCore.Qt.Checked:
            self.RHESSI_12_graph.setVisible(True)
        else:
            self.RHESSI_12_graph.setVisible(False)

    def RHESSI_3_clickBox(self, state):
        if state == QtCore.Qt.Checked:
            self.RHESSI_3_graph.setVisible(True)
        else:
            self.RHESSI_3_graph.setVisible(False)

    def RHESSI_6_clickBox(self, state):
        if state == QtCore.Qt.Checked:
            self.RHESSI_6_graph.setVisible(True)
        else:
            self.RHESSI_6_graph.setVisible(False)

    def RHESSI_50_clickBox(self, state):
        if state == QtCore.Qt.Checked:
            self.RHESSI_50_graph.setVisible(True)
        else:
            self.RHESSI_50_graph.setVisible(False)

    def GOES_XRSA_clickBox(self, state):
        if state == QtCore.Qt.Checked:
            self.GOES_XRSA_graph.setVisible(True)
        else:
            self.GOES_XRSA_graph.setVisible(False)

    def RHESSI_25_clickBox(self, state):
        if state == QtCore.Qt.Checked:
            self.RHESSI_25_graph.setVisible(True)
        else:
            self.RHESSI_25_graph.setVisible(False)

    def RHESSI_100_clickBox(self, state):
        if state == QtCore.Qt.Checked:
            self.RHESSI_100_graph.setVisible(True)
        else:
            self.RHESSI_100_graph.setVisible(False)

    def RHESSI_300_clickBox(self, state):
        if state == QtCore.Qt.Checked:
            self.RHESSI_300_graph.setVisible(True)
        else:
            self.RHESSI_300_graph.setVisible(False)

    def RHESSI_800_clickBox(self, state):
        if state == QtCore.Qt.Checked:
            self.RHESSI_800_graph.setVisible(True)
        else:
            self.RHESSI_800_graph.setVisible(False)

    def RHESSI_7000_clickBox(self, state):
        if state == QtCore.Qt.Checked:
            self.RHESSI_7000_graph.setVisible(True)
        else:
            self.RHESSI_7000_graph.setVisible(False)

    def show(self):
        self.main_window.show()

    def play_pause_func(self, play_pause_button):
        global playing
        if playing == True:
            playing = False
            print("stopped")
        else:
            playing = True
            print("start")

    def one_forward_func(self, one_forward_button):
        self.i += 1
        self.update_image()
        self.update_time_string()

    def one_backward_func(self, one_backward_button):
        self.i -= 1
        self.update_image()
        self.update_time_string()

    def speed_slider_func(self):
        self.frame_timer = 1000 - self.speed_slider.value()

    def movie_controls(self):
        play_pause_button = self.ui.play_pause_button
        play_pause_button_tab2 = self.ui.play_pause_button_tab2
        play_pause_button_tab3 = self.ui.play_pause_button_tab3
        one_forward_button = self.ui.one_forward_button
        one_forward_button_tab2 = self.ui.one_forward_button_tab2
        one_forward_button_tab3 = self.ui.one_forward_button_tab3
        one_backward_button = self.ui.one_backward_button
        one_backward_button_tab2 = self.ui.one_backward_button_tab2
        one_backward_button_tab3 = self.ui.one_backward_button_tab3
        self.change_hmi_view_button = self.ui.change_HMI_view_button
        self.change_aia_view_button = self.ui.change_AIA_view_button
        self.calculate_ribbon_area_button = self.ui.calculate_ribbon_area_button
        self.calculate_ribbon_area_button2 = self.ui.calculate_ribbon_area_button2
        self.calculate_ribbon_area_button3 = self.ui.calculate_ribbon_area_button3
        self.return_to_origin_button = self.ui.return_to_origin_button
        self.return_to_origin_button2 = self.ui.return_to_origin_button2
        self.return_to_origin_button3 = self.ui.return_to_origin_button3
        self.speed_slider = self.ui.speed_slider
        self.speed_slider.setMaximum(999)
        self.speed_slider.setMinimum(1)
        self.speed_slider.setValue(980)
        self.speed_slider_tab2 = self.ui.speed_slider_tab2
        self.speed_slider_tab2.setMaximum(999)
        self.speed_slider_tab2.setMinimum(1)
        self.speed_slider_tab2.setValue(980)
        self.speed_slider_tab3 = self.ui.speed_slider_tab3
        self.speed_slider_tab3.setMaximum(999)
        self.speed_slider_tab3.setMinimum(1)
        self.speed_slider_tab3.setValue(980)

        play_pause_button.clicked.connect(lambda: self.play_pause_func(play_pause_button))
        one_forward_button.clicked.connect(lambda: self.one_forward_func(one_forward_button))
        one_backward_button.clicked.connect(lambda: self.one_backward_func(one_backward_button))
        play_pause_button_tab2.clicked.connect(lambda: self.play_pause_func(play_pause_button_tab2))
        one_forward_button_tab2.clicked.connect(lambda: self.one_forward_func(one_forward_button_tab2))
        one_backward_button_tab2.clicked.connect(lambda: self.one_backward_func(one_backward_button_tab2))
        play_pause_button_tab3.clicked.connect(lambda: self.play_pause_func(play_pause_button_tab3))
        one_forward_button_tab3.clicked.connect(lambda: self.one_forward_func(one_forward_button_tab3))
        one_backward_button_tab3.clicked.connect(lambda: self.one_backward_func(one_backward_button_tab3))
        self.change_hmi_view_button.clicked.connect(lambda: self.change_hmi_view_func(self.change_hmi_view_button))
        self.change_aia_view_button.clicked.connect(lambda: self.change_aia_view_func(self.change_aia_view_button))
        self.speed_slider.valueChanged.connect(self.speed_slider_func)
        self.speed_slider_tab2.valueChanged.connect(self.speed_slider_func)
        self.speed_slider_tab3.valueChanged.connect(self.speed_slider_func)
        self.calculate_ribbon_area_button.clicked.connect(
            lambda: self.calculate_ribbon_area_func(self.calculate_ribbon_area_button))
        self.calculate_ribbon_area_button2.clicked.connect(
            lambda: self.calculate_ribbon_area_func(self.calculate_ribbon_area_button2))
        self.calculate_ribbon_area_button3.clicked.connect(
            lambda: self.calculate_ribbon_area_func(self.calculate_ribbon_area_button3))


    def calculate_ribbon_area_func(self, button):
        self.calculate_HMI_ribbons_centroid()
        data = all_wavelength_dict[1600]['full_raw_image_data'][self.i]
        header = all_wavelength_dict[1600]['header_list'][self.i]
        data = np.rot90(data)
        data = np.flipud(data)
        map_1600 = sunpy.map.Map(data, header)
        x, y = self.full_HMI_image.get_hmi_center()
        map_1600_cea, _ = HMI_to_CEA_func(map_1600, 2048, 4096)

        HMI_positive_array = copy.deepcopy(self.full_HMI_image.hmi_array)
        top_right = SkyCoord((x + 10) * u.deg, (y + 5) * u.deg,
                             frame=map_1600_cea.coordinate_frame)
        bottom_left = SkyCoord((x - 10) * u.deg,
                               (y - 5) * u.deg, frame=map_1600_cea.coordinate_frame)

        submap_1600 = (map_1600_cea.submap(bottom_left, top_right=top_right))

        negative_ribbon_area, positive_ribbon_area = self.full_HMI_image.calculate_ribbon_area(HMI_positive_array,
                                                                                               submap_1600,
                                                                                               int(self.cropped_1600_image_widget.max_clip_val))
        self.ui.positive_ribbon_area_value.setText(str(int(positive_ribbon_area)) + 'Mega-meters²')
        self.ui.negative_ribbon_area_value.setText(str(int(negative_ribbon_area)) + 'Mega-meters²')


    def return_to_origin_func(self):
        pass

    def change_aia_view_func(self, button):
        bruh = 300
        x, y = self.full_HMI_image.get_HMI_centroid_coords()
        x = int(x)
        y = int(y)
        for widget in self.full_widget_list:
            widget.draw_rectangle(x, y)
        self.initialize_aia_centering_label(x, y)

        for widget in self.cropped_widget_list:
            widget.x = x
            widget.y = y
            cropped_image_list = []

            if widget == self.cropped_aia_map_widget or widget == self.cropped_rhessi_map_widget:
                for array in all_wavelength_dict[widget.wavelength]['full_raw_image_data']:
                    cropped_image_data = array[x - bruh:x + bruh, y - bruh:y + bruh]
                    cropped_image_list.append(cropped_image_data)
                colormap = ct.aia_color_table(self.cropped_aia_map_widget.wavelength * u.angstrom)
                cropped_image_list = asinh_stretch_images(colormap, 0,
                                                          all_wavelength_dict[self.cropped_aia_map_widget.wavelength][
                                                              'max_clip_val'],
                                                          cropped_image_list)
                self.cropped_aia_map_widget.image_list = cropped_image_list

            else:
                widget.create_contours_list(x, y)

            widget.update_axis_dict(x, y)

            widget.change_cropped_aia_map_axis()

            widget.draw_cropped_limb()
            if widget.limb_state == True:
                widget.hide_limb()
                widget.show_limb_func()

            if widget.hmi_center_hide_or_show:
                print("widget shown")
                widget.hide_hmi_center()
                widget.show_hmi_center()

    def change_hmi_view_func(self, button):
        x, y = self.full_HMI_image.change_hmi_view_func()

    def update_time_string(self):
        for widget in self.widget_list:
            widget.time_label.setText(widget.time_string_list[self.i])

    def update_image(self):
        self.initialize_wavelength_and_state_values()

        if self.i > self.length_of_list:
            self.i = 0

        if self.i < 0:
            self.i = self.length_of_list

        for widget in self.widget_list:
            widget.update_aia_image(self.i)

        if self.cropped_aia_map_widget.show_center:
            self.cropped_aia_map_widget.show_center_mark(self.i)

        widget_contour_list_reference_dict = {
            self.cropped_aia_map_widget: [self.cropped_1600_image_widget.contour_display_list,
                                          self.cropped_helioprojective_HMI_image_widget.polarity_display_list],
            self.cropped_1600_image_widget: [self.cropped_1600_image_widget.contour_display_list1,
                                             self.cropped_helioprojective_HMI_image_widget.polarity_display_list1],
            self.cropped_helioprojective_HMI_image_widget: [self.cropped_1600_image_widget.contour_display_list2,
                                                            self.cropped_helioprojective_HMI_image_widget.polarity_display_list2],
            self.cropped_rhessi_map_widget: [self.cropped_1600_image_widget.contour_display_list3,
                                             self.cropped_helioprojective_HMI_image_widget.polarity_display_list3]
            }

        for widget in self.cropped_widget_list:
            if widget.contour_hide_or_show:
                widget.draw_1600_contours_on_widget(widget_contour_list_reference_dict[widget][0][self.i])
            if widget.polarities_hide_or_show:
                widget.draw_HMI_polarities_on_widget(
                    widget_contour_list_reference_dict[widget][1][self.i])

        self.real_time_update_centering_labels(self.i)
        self.lr.setValue(self.main_date_list[self.i])

    def run_widgets_func(self):
        self.updateTime = ptime.time()

        self.fps = 0

        def updateData():
            if playing:
                self.i = self.i + 1
                self.update_image()
                self.update_time_string()

            QTimer.singleShot(self.frame_timer, updateData)

            # getting current time
            now = ptime.time()

            # temporary fps
            fps2 = 1.0 / (now - self.updateTime)

            # updating the time
            self.updateTime = now

            # setting original fps value
            self.fps = self.fps * 0.9 + fps2 * 0.1

        updateData()

    def initialize_time_list(self):
        list1 = all_wavelength_dict[171]['date_list']
        self.main_date_list = []
        self.alt_date_list = []
        for date1 in list1:
            a = (date1 - datetime(1969, 12, 31, 17, 0, 0)).total_seconds()
            self.main_date_list.append(a)

    def goes_euv_plot(self):

        night_list = self.get_night_saa_time_range(special_beginning_time_str, special_end_time_str,
                                                   'Night_time_data.csv')
        SAA_list = self.get_night_saa_time_range(special_beginning_time_str, special_end_time_str, 'SAA_time_data.csv')

        pw = self.ui.graphs_widget

        self.GOES_xrsb = self.ui.GOES_XRSB
        self.GOES_xrsa = self.ui.GOES_XRSA
        self.RHESSI_12 = self.ui.RHESSI_12
        self.RHESSI_6 = self.ui.RHESSI_6
        self.RHESSI_50 = self.ui.RHESSI_50
        self.RHESSI_3 = self.ui.RHESSI_3
        self.RHESSI_25 = self.ui.RHESSI_25
        self.RHESSI_100 = self.ui.RHESSI_100
        self.RHESSI_300 = self.ui.RHESSI_300
        self.RHESSI_800 = self.ui.RHESSI_800
        self.RHESSI_7000 = self.ui.RHESSI_7000

        self.EUV_graph_list = [self.ui.EUV_flux_94, self.ui.EUV_flux_171,  self.ui.EUV_flux_211,
                               self.ui.EUV_flux_304, self.ui.EUV_flux_1600]
        self.normalized_EUV_graph_list = [self.ui.normalized_EUV_flux_94, self.ui.normalized_EUV_flux_171,
                                          self.ui.normalized_EUV_flux_211,
                                          self.ui.normalized_EUV_flux_304, self.ui.normalized_EUV_flux_1600]

        self.GOES_xrsb.stateChanged.connect(self.GOES_XRSB_clickBox)
        self.GOES_xrsa.stateChanged.connect(self.GOES_XRSA_clickBox)
        self.RHESSI_50.stateChanged.connect(self.RHESSI_50_clickBox)
        self.RHESSI_3.stateChanged.connect(self.RHESSI_3_clickBox)
        self.RHESSI_12.stateChanged.connect(self.RHESSI_12_clickBox)

        self.RHESSI_25.stateChanged.connect(self.RHESSI_25_clickBox)
        self.RHESSI_50.stateChanged.connect(self.RHESSI_50_clickBox)
        self.RHESSI_100.stateChanged.connect(self.RHESSI_100_clickBox)
        self.RHESSI_300.stateChanged.connect(self.RHESSI_300_clickBox)
        self.RHESSI_800.stateChanged.connect(self.RHESSI_800_clickBox)
        self.RHESSI_7000.stateChanged.connect(self.RHESSI_7000_clickBox)
        self.RHESSI_6.stateChanged.connect(self.RHESSI_6_clickBox)
        self.RHESSI_7000.stateChanged.connect(self.RHESSI_7000_clickBox)

        goes_pickle = 'Patched_GOES_Data/' + event_beginning_time1[:4] + '_GOES_DATA.pickle'
        RHESSI_pickle = event_beginning_time1[:4] + '_final_RHESSI_flux.pickle'
        dataframe = pd.read_pickle(goes_pickle)[self.beginning_time:self.end_time]

        time = list(dataframe.index)
        xrsa = dataframe['xrsa']
        xrsb = list(dataframe['xrsb'])
        self.new_time = []

        RHESSI_dataframe = pd.read_pickle(RHESSI_pickle)[self.beginning_time:self.end_time]
        RHESSI_wavelength_list = ['3 - 6 keV', '6 - 12 keV', '12 - 25 keV', '25 - 50 keV', '50 - 100 keV',
                                  '100 - 300 keV', '300 - 800 keV', '800 - 7000 keV', '7000 - 20000 keV']

        RHESSI_time_list = list(RHESSI_dataframe.index)
        RHESSI_timeseries_dict = {}

        for wavelength in RHESSI_wavelength_list:
            RHESSI_timeseries_dict[wavelength] = list(RHESSI_dataframe[wavelength])

        RHESSI_final_time_list = []

        for date1 in time:
            a = (date1 - datetime(1969, 12, 31, 17, 0, 0)).total_seconds()
            self.new_time.append(a)

        for date1 in RHESSI_time_list:
            a = (date1 - datetime(1969, 12, 31, 17, 0, 0)).total_seconds()
            RHESSI_final_time_list.append(a)

        pi_1 = pw.plotItem
        self.vb_1 = pi_1.getViewBox()  # use original viewbox
        self.vb_1.setBorder('#777')
        self.vb_2 = pg.ViewBox()  # prepare additional viewbox 2
        self.vb_3 = pg.ViewBox()
        pi_1.scene().addItem(self.vb_2)
        pi_1.scene().addItem(self.vb_3)

        # we are putting all y-axes on the right, so we do not need this one.
        ax_1 = pi_1.getAxis('left')  # use original right axis of pi_1
        ax_1.setLabel('Watts/m⁻²', color='#e33')
        ax_1.enableAutoSIPrefix(enable=False)

        ax_1.setLogMode(True)
        # ax_1.show() # the right side axis is hidden by default

        self.ax_2 = pg.AxisItem('right')  # prepare axis 2...
        self.ax_2.setLabel('EUV Flux (DN/sec)', color='#3e3')
        self.ax_2.enableAutoSIPrefix(enable=False)
        self.ax_2.setLogMode(False)
        self.ax_2.linkToView(self.vb_2)  # ...assign it to viewbox 2...
        pi_1.layout.addItem(self.ax_2, 2, 3)  # ...and add it to the right of the original axis

        ax_3 = pg.AxisItem('right')  # prepare axis 3...
        ax_3.setLabel('Total Counts', color='#33e')
        ax_3.setLogMode(True)
        ax_3.enableAutoSIPrefix(enable=False)
        ax_3.linkToView(self.vb_3)  # ...assign it to viewbox 3...
        pi_1.layout.addItem(ax_3, 2, 4)  # ...and add it to the right of the second axis

        # --- Handle view resizing ---
        def updateViews():
            ## view has resized; update auxiliary views to match
            self.vb_2.setGeometry(self.vb_1.sceneBoundingRect())
            self.vb_3.setGeometry(self.vb_1.sceneBoundingRect())
            ## update linked axes, may not be necessary anymore.
            ## (probably this should be handled in ViewBox.resizeEvent)
            # vb_2.linkedViewChanged(vb_1, vb_2.XAxis)
            # vb_3.linkedViewChanged(vb_1, vb_3.XAxis)

        updateViews()  # make sure all viewboxes are drawn in the same place...
        self.vb_1.sigResized.connect(updateViews)  # ... and get realigned after any resize.

        # --- Add plot data ---
        self.GOES_XRSA_graph = pg.PlotDataItem(self.new_time, xrsa, pen=pg.mkPen('b', width=1))
        self.GOES_XRSA_graph.setLogMode(False, True)  # y is linear
        self.GOES_xrsa.setChecked(True)
        self.vb_1.addItem(self.GOES_XRSA_graph)

        self.GOES_XRSB_graph = pg.PlotDataItem(self.new_time, xrsb, pen=pg.mkPen('r', width=1))
        self.GOES_XRSB_graph.setLogMode(False, True)  # y is linear
        self.GOES_xrsb.setChecked(True)
        self.vb_1.addItem(self.GOES_XRSB_graph)
        """
        self.EUV_flux_171_graph = pg.PlotDataItem(self.main_date_list, all_wavelength_dict[171]['intensity_list'],
                                                  pen=pg.mkPen('g', width=1))
        # initialize_EUV_plots(self.EUV_flux_171_graph)

        self.EUV_flux_171_graph.setLogMode(False, False)  # y is logarithmic
        self.EUV_flux_171.setChecked(True)
        self.vb_2.addItem(self.EUV_flux_171_graph)
        self.vb_2.setXLink(vb_1)
        """
        self.EUV_wavelength_list =  [94, 171, 211, 304, 1600]
        self.EUV_graph_dict = {}
        self.normalized_EUV_graph_dict = {}
        color_list = ['r', 'g', 'b', 'y', 'c', 'l']
        max_list = []
        min_list = []

        for i in range(0, len(self.EUV_graph_list)):
            intensity_list = all_wavelength_dict[self.EUV_wavelength_list[i]]['intensity_list']
            graph = Graph(date_list=all_wavelength_dict[self.EUV_wavelength_list[i]]['date_list'],
                          intensity_list=intensity_list, vb_1=self.vb_1, vb_2=self.vb_2,
                          clickbox=self.EUV_graph_list[i], color=color_list[i])

            self.EUV_graph_dict[self.EUV_wavelength_list[i]] = graph

            normalized_graph = Graph(date_list=all_wavelength_dict[self.EUV_wavelength_list[i]]['date_list'],
                                     intensity_list=all_wavelength_dict[self.EUV_wavelength_list[i]][
                                         'normalized_intensity_list'], vb_1=self.vb_1,
                                     vb_2=self.vb_2, clickbox=self.normalized_EUV_graph_list[i], color=color_list[i])
            # normalized_graph.graph.setVisible(False)
            self.normalized_EUV_graph_dict[self.EUV_wavelength_list[i]] = normalized_graph
            max_list.append(max(intensity_list))
            min_list.append(min(intensity_list))

        self.min_intensity_val = min(min_list)
        self.max_intensity_val = max(max_list)


        self.RHESSI_3_graph = pg.PlotDataItem(RHESSI_final_time_list, RHESSI_timeseries_dict['3 - 6 keV'],
                                              pen=pg.mkPen('b', width=1))
        self.RHESSI_3_graph.setLogMode(False, True)
        self.vb_3.addItem(self.RHESSI_3_graph)
        self.RHESSI_3.setChecked(True)
        self.RHESSI_6_graph = pg.PlotDataItem(RHESSI_final_time_list, RHESSI_timeseries_dict['6 - 12 keV'],
                                              pen=pg.mkPen('c', width=1))
        self.RHESSI_6_graph.setLogMode(False, True)
        self.vb_3.addItem(self.RHESSI_6_graph)
        self.RHESSI_6.setChecked(True)
        self.RHESSI_12_graph = pg.PlotDataItem(RHESSI_final_time_list, RHESSI_timeseries_dict['12 - 25 keV'],
                                               pen=pg.mkPen('m', width=1))
        self.RHESSI_12_graph.setLogMode(False, True)
        self.RHESSI_12.setChecked(True)
        self.vb_3.addItem(self.RHESSI_12_graph)
        self.RHESSI_25_graph = pg.PlotDataItem(RHESSI_final_time_list, RHESSI_timeseries_dict['25 - 50 keV'],
                                               pen=pg.mkPen('y', width=1))
        self.RHESSI_25_graph.setLogMode(False, True)
        self.RHESSI_25.setChecked(True)
        self.vb_3.addItem(self.RHESSI_25_graph)
        self.RHESSI_50_graph = pg.PlotDataItem(RHESSI_final_time_list, RHESSI_timeseries_dict['50 - 100 keV'],
                                               pen=pg.mkPen('k', width=1))
        self.RHESSI_50_graph.setLogMode(False, True)
        self.RHESSI_50.setChecked(True)
        self.vb_3.addItem(self.RHESSI_50_graph)
        self.RHESSI_100_graph = pg.PlotDataItem(RHESSI_final_time_list, RHESSI_timeseries_dict['100 - 300 keV'],
                                                pen=pg.mkPen('w', width=1))
        self.RHESSI_100_graph.setLogMode(False, True)
        self.RHESSI_100.setChecked(True)
        self.vb_3.addItem(self.RHESSI_100_graph)
        self.RHESSI_300_graph = pg.PlotDataItem(RHESSI_final_time_list, RHESSI_timeseries_dict['300 - 800 keV'],
                                                pen=pg.mkPen('r', width=1))
        self.RHESSI_300_graph.setLogMode(False, True)
        self.RHESSI_300.setChecked(True)
        self.vb_3.addItem(self.RHESSI_300_graph)
        self.RHESSI_800_graph = pg.PlotDataItem(RHESSI_final_time_list, RHESSI_timeseries_dict['800 - 7000 keV'],
                                                pen=pg.mkPen('g', width=1))
        self.RHESSI_800_graph.setLogMode(False, True)
        self.RHESSI_800.setChecked(True)
        self.vb_3.addItem(self.RHESSI_800_graph)
        self.RHESSI_7000_graph = pg.PlotDataItem(RHESSI_final_time_list, RHESSI_timeseries_dict['7000 - 20000 keV'],
                                                 pen=pg.mkPen('b', width=1))

        self.RHESSI_7000_graph.setLogMode(False, True)
        self.RHESSI_7000.setChecked(True)
        self.vb_3.addItem(self.RHESSI_7000_graph)

        self.vb_3.setXLink(self.vb_1)

        self.lr = pg.InfiniteLine()
        self.lr.setMovable(True)
        self.lr.setValue(self.main_date_list[0])
        self.lr.setZValue(1)

        pi_1.addItem(self.lr)

        for time_range in night_list:
            linear_region = pg.LinearRegionItem(time_range, pen=(0, 0, 0, 0), brush=(128, 255, 0, 60))
            linear_region.setMovable(False)
            linear_region.setZValue(1)
            pi_1.addItem(linear_region)

        for time_range in SAA_list:
            linear_region = pg.LinearRegionItem(time_range, pen=(0, 0, 0, 0), brush=(255, 0, 0, 60))
            linear_region.setMovable(False)
            linear_region.setZValue(1)
            pi_1.addItem(linear_region)

        pw.setXRange(self.main_date_list[0], self.main_date_list[len(self.main_date_list) - 1])

        def handle_sig_dragged(obj):
            assert obj is self.lr
            value = take_closest(self.main_date_list, obj.value())
            index = self.main_date_list.index(value)
            self.i = index
            self.update_image()
            self.update_time_string()

        self.lr.sigDragged.connect(handle_sig_dragged)

    def full_1600_image_widget_func(self):
        self.cropped_aia_1600_list = all_wavelength_dict[1600]['cropped_image_data']
        self.full_aia_1600_list = all_wavelength_dict[1600]['image_data']
        self.date_list_1600 = all_wavelength_dict[1600]['time_string_list']


if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec())

