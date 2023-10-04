from units_conversion_functions import heliographic_to_helioprojective
from astropy.coordinates import SkyCoord
from reproject import reproject_interp
import sunpy.map
from pffspy_util import car_to_cea
import copy
import astropy.wcs
from astropy.wcs import WCS
import numpy as np
import astropy.units as u


def unlock_heliographic_coords(x, y, m_hmi_cea):
    heliographic_coord = WCS(m_hmi_cea.meta).pixel_to_world(x, y)
    print(heliographic_coord)
    NS_lat = float(str(heliographic_coord.lat).replace('.', '').replace('d', '.').replace('m', '').replace('s', ''))
    EW_lon = float(str(heliographic_coord.lon).replace('.', '').replace('d', '.').replace('m', '').replace('s', ''))
    print(NS_lat, EW_lon)
    print("coords")
    return NS_lat, EW_lon


def HMI_weighted_centering_func(array1, string):
    new_array = array1
    print(new_array.shape)
    print(new_array)
    if string == "greater":
        print("greater")
        new_array[new_array < 50] = 0
        print(new_array)
    if string == "less":
        print("less")

        new_array = new_array * -1

        new_array[new_array < 50] = 0

    x_val = range(0, new_array.shape[0])
    y_val = range(0, new_array.shape[1])

    (X, Y) = np.meshgrid(x_val, y_val)

    x_coord = (X * new_array).sum() / new_array.sum().astype("float")
    y_coord = (Y * new_array).sum() / new_array.sum().astype("float")

    return [x_coord, y_coord]

def centeroidnp(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])

    return sum_x / length, sum_y / length


def get_HMI_data(x,y):
    sun_center_x_list_placeholder = []
    sun_center_y_list_placeholder = []
    pixel_arcsec_x_list_placeholder = []
    pixel_arcsec_y_list_placeholder = []


        # return 1,2

    m_hmi = sunpy.map.Map('/Users/jamisenma/sunpy/data/hmi_m_45s_2011_02_14_02_00_00_tai_magnetogram.fits')

    print("cooridnate frame")

    print(m_hmi.coordinate_frame)

    shape_out = [2048, 4096]
    frame_out = SkyCoord(0, 0, unit=u.deg, rsun=m_hmi.coordinate_frame.rsun, frame="heliographic_stonyhurst",
                         obstime=m_hmi.date)
    header = sunpy.map.make_fitswcs_header(
        shape_out,
        frame_out,
        scale=[180 / shape_out[0], 360 / shape_out[1]] * u.deg / u.pix,
        projection_code="CAR",
    )
    print(shape_out[0], shape_out[1], u.deg, u.pix)
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
    m_hmi_cea_copy = copy.deepcopy(m_hmi_cea)
    print("header")
    print(m_hmi_cea.meta)


    full_data = (m_hmi_cea.data)
    full_data2 = copy.deepcopy(full_data)


    print("meta hmi")
    print(m_hmi_cea.meta['crpix1'])
    sun_center_x_list_placeholder.append(m_hmi_cea.meta['crpix1'])
    sun_center_y_list_placeholder.append(m_hmi_cea.meta['crpix2'])
    pixel_arcsec_x_list_placeholder.append(m_hmi_cea.meta['cdelt1'])
    pixel_arcsec_y_list_placeholder.append(m_hmi_cea.meta['cdelt2'])



    #x = 2000
    #y = 2048 - 1369



    side_length = 300


    positive_array = full_data[y - side_length:y + side_length,
                     x - side_length:x + side_length]

    negative_array = copy.deepcopy(positive_array)

    x_positive, y_positive = HMI_weighted_centering_func(positive_array, 'greater')

    x_negative, y_negative = HMI_weighted_centering_func(negative_array, 'less')
    # print(x_positive,y_positive)

    x_positive = x_positive + x - 300
    x_negative = x_negative + x - 300
    y_negative = y_negative + y - 300
    y_positive = y_positive + y - 300

    x_negative_heliographic, y_negative_heliographic = unlock_heliographic_coords(x_negative, y_negative, m_hmi_cea)
    x_positive_heliographic, y_positive_heliographic = unlock_heliographic_coords(x_positive, y_positive, m_hmi_cea)
    hmi_heliographic_x = (x_positive_heliographic + x_negative_heliographic) / 2
    hmi_heliographic_y = (y_positive_heliographic + y_negative_heliographic) / 2
    x_positive_helioprojective, y_positive_helioprojective = heliographic_to_helioprojective(hmi_heliographic_x,
                                                                                             hmi_heliographic_y,
                                                                                             '2011-02-15')

    hmi_center_x_pixel = (x_negative + x_positive) / 2
    hmi_center_y_pixel = (y_negative + y_positive) / 2

    return m_hmi_cea, x_negative, x_positive, y_negative, y_positive, hmi_center_x_pixel, hmi_center_y_pixel, full_data2, m_hmi_cea_copy, out_wcs

def get_pixel_coords_from_AR_center(region, date1):
    AR_number = '1' + str(region['AR'])
    location = str(region['Location'])
    hale_class = str(region['Hale_Class'])

    print(location)
    if 'N' in location:

        NS = 1
    else:
        NS = -1

    if 'W' in location:
        location = location.split('W')
        EW = 1
    else:
        EW = -1
        location = location.split('E')

    print(location)
    NS = NS * int(location[0][1:])
    EW = EW * int(location[1])

    x_helioprojective, y_helioprojective = heliographic_to_helioprojective(NS, EW, date1)
    x_pixel = x_helioprojective_to_pixel(x_helioprojective, 0,  pixel_arcsec_x_list_placeholder[0], sun_center_x_list_placeholder[0]) / 8
    y_pixel = y_helioprojective_to_pixel(y_helioprojective, 0, pixel_arcsec_y_list_placeholder[0], sun_center_y_list_placeholder[0]) / 8
    return x_pixel, y_pixel, AR_number, hale_class