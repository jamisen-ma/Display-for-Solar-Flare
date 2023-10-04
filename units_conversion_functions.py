import astropy.units as u
from astropy.coordinates import SkyCoord
import re
from sunpy.coordinates import frames

def heliographic_to_helioprojective(x, y, time):
    c = SkyCoord(y * u.deg, x * u.deg, frame=frames.HeliographicStonyhurst, obstime=time, observer="earth")

    helioprojective_coord = (c.transform_to(frames.Helioprojective))

    print(helioprojective_coord)

    NS_lat = float(re.sub(r'[a-z]+', '', str(helioprojective_coord.Tx)[:12], re.I))
    EW_lon = float(re.sub(r'[a-z]+', '', str(helioprojective_coord.Ty)[:12], re.I))

    return NS_lat, EW_lon
    # return 1,2


def x_helioprojective_to_pixel(value, additional_value, pixel_arcsec_x_list, sun_center_x_list):
    return ((value + additional_value) / pixel_arcsec_x_list[0]) + sun_center_x_list[0]


def y_helioprojective_to_pixel(value, additional_value, pixel_arcsec_y_list, sun_center_y_list):
    return ((-value + additional_value) / pixel_arcsec_y_list[0]) + sun_center_y_list[0]


def x_pixel_to_helioprojective(value, sun_center_x_list, pixel_arcsec_x_list):
    return (value - sun_center_x_list[0]) * pixel_arcsec_x_list[0]


def y_pixel_to_helioprojective(value, sun_center_y_list, pixel_arcsec_y_list):
    return -(value - sun_center_y_list[0]) * pixel_arcsec_y_list[0]


def rhessi_pixel_to_helioprojective_x(header_center, value, multiplier):
    return (value - header_center) * multiplier


def rhessi_pixel_to_helioprojective_y(header_center, value, multiplier):
    return -(value - header_center) * multiplier


def rhessi_x_helioprojective_to_pixel(header_center, value, multiplier, additional_value):
    return ((value + additional_value) / multiplier) + header_center


def rhessi_y_helioprojective_to_pixel(header_center, value, multiplier, additional_value):
    return ((-value + additional_value) / multiplier) + header_center


def helioprojective_to_heliographic_simplified(x, y, time):
    c = SkyCoord(x * u.arcsec, y * u.arcsec, frame=frames.Helioprojective, obstime=time, observer="earth")

    heliographic_coord = (c.transform_to(frames.HeliographicStonyhurst))

    print(heliographic_coord)

    # North and South is latitude
    # West and East is longitude

    NS_lat = int(str(heliographic_coord.lat).split('d')[0])
    EW_lon = int(str(heliographic_coord.lon).split('d')[0])
    print(NS_lat, EW_lon)
    return NS_lat, EW_lon


def helioprojective_to_heliographic(x, y, time):
    NS_lat, EW_lon = helioprojective_to_heliographic_simplified(x, y, time)
    if NS_lat > 0:
        first_part = "N" + str(NS_lat)
    else:
        first_part = "S" + str(abs(NS_lat))

    if EW_lon > 0:
        second_part = "W" + str(EW_lon)
    else:
        second_part = "E" + str(abs(EW_lon))
    return first_part + second_part