from __future__ import annotations

import re
from typing import List, Tuple

import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames

from logging_setup import get_logger


logger = get_logger(__name__)


def heliographic_to_helioprojective(x: float, y: float, time: str) -> Tuple[float, float]:
    c = SkyCoord(y * u.deg, x * u.deg, frame=frames.HeliographicStonyhurst, obstime=time, observer="earth")
    helioprojective_coord = c.transform_to(frames.Helioprojective)
    NS_lat = float(re.sub(r"[a-z]+", "", str(helioprojective_coord.Tx)[:12], re.I))
    EW_lon = float(re.sub(r"[a-z]+", "", str(helioprojective_coord.Ty)[:12], re.I))
    return NS_lat, EW_lon


def x_helioprojective_to_pixel(value: float, additional_value: float, pixel_arcsec_x_list: List[float], sun_center_x_list: List[float]) -> float:
    return ((value + additional_value) / pixel_arcsec_x_list[0]) + sun_center_x_list[0]


def y_helioprojective_to_pixel(value: float, additional_value: float, pixel_arcsec_y_list: List[float], sun_center_y_list: List[float]) -> float:
    return ((-value + additional_value) / pixel_arcsec_y_list[0]) + sun_center_y_list[0]


def x_pixel_to_helioprojective(value: float, sun_center_x_list: List[float], pixel_arcsec_x_list: List[float]) -> float:
    return (value - sun_center_x_list[0]) * pixel_arcsec_x_list[0]


def y_pixel_to_helioprojective(value: float, sun_center_y_list: List[float], pixel_arcsec_y_list: List[float]) -> float:
    return -(value - sun_center_y_list[0]) * pixel_arcsec_y_list[0]


def rhessi_pixel_to_helioprojective_x(header_center: float, value: float, multiplier: float) -> float:
    return (value - header_center) * multiplier


def rhessi_pixel_to_helioprojective_y(header_center: float, value: float, multiplier: float) -> float:
    return -(value - header_center) * multiplier


def rhessi_x_helioprojective_to_pixel(header_center: float, value: float, multiplier: float, additional_value: float) -> float:
    return ((value + additional_value) / multiplier) + header_center


def rhessi_y_helioprojective_to_pixel(header_center: float, value: float, multiplier: float, additional_value: float) -> float:
    return ((-value + additional_value) / multiplier) + header_center


def helioprojective_to_heliographic_simplified(x: float, y: float, time: str) -> Tuple[int, int]:
    c = SkyCoord(x * u.arcsec, y * u.arcsec, frame=frames.Helioprojective, obstime=time, observer="earth")
    heliographic_coord = c.transform_to(frames.HeliographicStonyhurst)
    NS_lat = int(str(heliographic_coord.lat).split('d')[0])
    EW_lon = int(str(heliographic_coord.lon).split('d')[0])
    return NS_lat, EW_lon


def helioprojective_to_heliographic(x: float, y: float, time: str) -> str:
    NS_lat, EW_lon = helioprojective_to_heliographic_simplified(x, y, time)
    first_part = ("N" if NS_lat > 0 else "S") + str(abs(NS_lat))
    second_part = ("W" if EW_lon > 0 else "E") + str(abs(EW_lon))
    return first_part + second_part