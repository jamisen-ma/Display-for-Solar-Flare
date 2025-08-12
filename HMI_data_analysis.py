from __future__ import annotations

from typing import Tuple

import copy
import numpy as np
import astropy.units as u
import astropy.wcs
import sunpy.map
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from reproject import reproject_interp

from pffspy_util import car_to_cea
from units_conversion_functions import (
    heliographic_to_helioprojective,
    x_helioprojective_to_pixel,
    y_helioprojective_to_pixel,
)
from constants import (
    sun_center_x_list_placeholder,
    sun_center_y_list_placeholder,
    pixel_arcsec_x_list_placeholder,
    pixel_arcsec_y_list_placeholder,
)
from logging_setup import get_logger


logger = get_logger(__name__)


def unlock_heliographic_coords(x: float, y: float, m_hmi_cea: sunpy.map.Map) -> Tuple[float, float]:
    heliographic_coord = WCS(m_hmi_cea.meta).pixel_to_world(x, y)
    NS_lat = float(str(heliographic_coord.lat).replace('.', '').replace('d', '.').replace('m', '').replace('s', ''))
    EW_lon = float(str(heliographic_coord.lon).replace('.', '').replace('d', '.').replace('m', '').replace('s', ''))
    return NS_lat, EW_lon


def HMI_weighted_centering_func(array1: np.ndarray, string: str):
    new_array = array1 * (string == "less") * -1
    new_array[new_array < 50] = 0
    X, Y = np.meshgrid(range(new_array.shape[0]), range(new_array.shape[1]))
    x_coord = (X * new_array).sum() / new_array.sum().astype("float")
    y_coord = (Y * new_array).sum() / new_array.sum().astype("float")
    return [x_coord, y_coord]


def get_HMI_data(x: int, y: int, fits_path: str = '/Users/jamisenma/sunpy/data/hmi_m_45s_2011_02_14_02_00_00_tai_magnetogram.fits'):
    m_hmi = sunpy.map.Map(fits_path)
    shape_out = [2048, 4096]
    frame_out = SkyCoord(0, 0, unit=u.deg, rsun=m_hmi.coordinate_frame.rsun, frame="heliographic_stonyhurst", obstime=m_hmi.date)
    header = sunpy.map.make_fitswcs_header(shape_out, frame_out, scale=[180 / shape_out[0], 360 / shape_out[1]] * u.deg / u.pix, projection_code="CAR")
    out_wcs = WCS(header)
    array, _ = reproject_interp(m_hmi, out_wcs, shape_out=shape_out)
    array = np.where(np.isnan(array), 0, array)
    m_hmi_cea = car_to_cea(sunpy.map.Map(array, header))
    m_hmi_cea.meta.update({'TELESCOP': m_hmi.meta['TELESCOP'], 'CONTENT': 'Carrington Synoptic Chart Of Br Field', 'T_OBS': m_hmi_cea.meta.pop('DATE-OBS')})
    m_hmi_cea = sunpy.map.Map(m_hmi_cea.data, m_hmi_cea.meta)

    full_data = m_hmi_cea.data
    # Placeholder geometry from constants
    # Values from map meta are retained locally if needed

    side_length = 300
    positive_array = full_data[y - side_length:y + side_length, x - side_length:x + side_length]
    x_positive, y_positive = HMI_weighted_centering_func(positive_array, 'greater')
    x_negative, y_negative = HMI_weighted_centering_func(positive_array, 'less')
    x_positive += x - 300
    x_negative += x - 300
    y_negative += y - 300
    y_positive += y - 300

    x_neg_hg, y_neg_hg = unlock_heliographic_coords(x_negative, y_negative, m_hmi_cea)
    x_pos_hg, y_pos_hg = unlock_heliographic_coords(x_positive, y_positive, m_hmi_cea)
    hmi_hg_x = (x_pos_hg + x_neg_hg) / 2
    hmi_hg_y = (y_pos_hg + y_neg_hg) / 2
    _x_pos_hp, _y_pos_hp = heliographic_to_helioprojective(hmi_hg_x, hmi_hg_y, '2011-02-15')

    hmi_center_x_pixel = (x_negative + x_positive) / 2
    hmi_center_y_pixel = (y_negative + y_positive) / 2

    return m_hmi_cea, x_negative, x_positive, y_negative, y_positive, hmi_center_x_pixel, hmi_center_y_pixel, full_data, copy.deepcopy(m_hmi_cea), out_wcs


def get_pixel_coords_from_AR_center(region, date1: str):
    AR_number = '1' + str(region['AR'])
    location = str(region['Location'])
    hale_class = str(region['Hale_Class'])
    NS = 1 if 'N' in location else -1
    EW = 1 if 'W' in location else -1
    loc_parts = location.replace('N', '').replace('S', '').replace('W', '').replace('E', '').split()
    NS *= int(loc_parts[0])
    EW *= int(loc_parts[1])
    x_hp, y_hp = heliographic_to_helioprojective(NS, EW, date1)
    x_pixel = x_helioprojective_to_pixel(x_hp, 0, pixel_arcsec_x_list_placeholder, sun_center_x_list_placeholder) / 8
    y_pixel = y_helioprojective_to_pixel(y_hp, 0, pixel_arcsec_y_list_placeholder, sun_center_y_list_placeholder) / 8
    return x_pixel, y_pixel, AR_number, hale_class
