from __future__ import annotations

from typing import Iterable, List, Sequence, Tuple

import numpy as np
from astropy.visualization.mpl_normalize import ImageNormalize
from image_processing_functions import convert_to_3D_array, strided_rescale
from astropy.visualization import AsinhStretch
import cv2
from datetime import date, timedelta
from sunpy.visualization.colormaps import color_tables as ct
import astropy.units as u

from logging_setup import get_logger


logger = get_logger(__name__)


def creating_1600_image_list(
    list_1600: Sequence[np.ndarray], x: int, y: int, exposure_list: Sequence[float]
) -> Tuple[List[np.ndarray], List[np.ndarray], List[float], List[list]]:
    peak_x = x
    peak_y = 4096 - y
    crop_half_size = 300
    sdoaia = ct.aia_color_table(1600 * u.angstrom)
    cropped_images: List[np.ndarray] = []
    full_1600_image_list: List[np.ndarray] = []
    intensity_list: List[float] = []
    contour_list: List[list] = []

    for i in range(0, len(list_1600)):
        data = np.clip(list_1600[i], 0, 1500)
        cropped = data[peak_x - crop_half_size : peak_x + crop_half_size,
                       peak_y - crop_half_size : peak_y + crop_half_size]
        intensity_value = float(np.divide(cropped.sum(), exposure_list[i]))
        intensity_list.append(intensity_value)
        cropped = ImageNormalize(stretch=AsinhStretch(0.01)).__call__(cropped)
        cropped_color = (np.uint8(sdoaia(cropped) * 255))

        low_res = strided_rescale(data, 8)
        low_res = ImageNormalize(stretch=AsinhStretch(0.01)).__call__(low_res)
        low_res_color = (np.uint8(sdoaia(low_res) * 255))

        cropped_images.append(cropped_color)
        full_1600_image_list.append(low_res_color)

    final_1600_list: List[np.ndarray] = []
    base_image = cropped_images[0]
    zero_frame = np.array(convert_to_3D_array(base_image.astype(np.int16)))
    for k in range(1, len(cropped_images)):
        img = cropped_images[k]
        one_frame = np.array(convert_to_3D_array(img).astype(np.int16))
        im = np.absolute(one_frame - zero_frame).astype(np.uint8)
        imgray = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)

        _, thresh = cv2.threshold(imgray, 80, 255, cv2.THRESH_BINARY)
        contours, _ = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
        cnts = sorted(contours, key=cv2.contourArea, reverse=True)[:5]
        number_of_contours = 4
        for j in range(min(number_of_contours, len(cnts))):
            cv2.drawContours(img, cnts[j], -1, (0, 255, 0), 1)

        final_1600_list.append(img)
        contour_list.append(cnts)

    return final_1600_list, full_1600_image_list, intensity_list, contour_list


def create_difference_and_ratio_lists(image_data_list: Sequence[np.ndarray], sdoaia):
    run_diff_list: List[np.ndarray] = []
    base_diff_list: List[np.ndarray] = []
    run_ratio_list: List[np.ndarray] = []
    base_ratio_list: List[np.ndarray] = []
    zero_image = image_data_list[0]
    for i in range(1, len(image_data_list)):
        first_image = image_data_list[i]
        second_image = image_data_list[i - 1]

        run_diff_image = strided_rescale(np.subtract(first_image, second_image), 8)
        run_diff_image = np.clip(run_diff_image, 0, 1500)
        run_diff_image = ImageNormalize(stretch=AsinhStretch(0.01)).__call__(run_diff_image)
        run_diff_image = (np.uint8(sdoaia(run_diff_image) * 255))

        base_diff_image = strided_rescale(np.subtract(first_image, zero_image), 8)
        base_diff_image = np.clip(base_diff_image, 0, 1500)
        base_diff_image = ImageNormalize(stretch=AsinhStretch(0.01)).__call__(base_diff_image)
        base_diff_image = (np.uint8(sdoaia(base_diff_image) * 255))

        run_ratio_image = strided_rescale(np.divide(first_image, second_image), 8)
        run_ratio_image = np.clip(run_ratio_image, 1, 10)
        run_ratio_image = ImageNormalize(stretch=AsinhStretch(0.01)).__call__(run_ratio_image)
        run_ratio_image = (np.uint8(sdoaia(run_ratio_image) * 255))

        base_ratio_image = strided_rescale(np.divide(first_image, zero_image), 8)
        base_ratio_image = np.clip(base_ratio_image, 1, 10)
        base_ratio_image = ImageNormalize(stretch=AsinhStretch(0.01)).__call__(base_ratio_image)
        base_ratio_image = (np.uint8(sdoaia(base_ratio_image) * 255))

        run_diff_list.append(run_diff_image)
        base_diff_list.append(base_diff_image)
        run_ratio_list.append(run_ratio_image)
        base_ratio_list.append(base_ratio_image)

    return run_diff_list, base_diff_list, run_ratio_list, base_ratio_list


def create_date_list(sdate: str, edate: str) -> List[str]:
    """Create repeated date list in YYYY-MM-DD across a start-end inclusive range.

    Returns each date twice to match expected consumer behavior.
    """
    start = date(int(sdate[:4]), int(sdate[5:7]), int(sdate[8:10]))
    end = date(int(edate[:4]), int(edate[5:7]), int(edate[8:10]))

    delta = end - start
    date_list: List[str] = []

    for i in range(delta.days):
        day = start + timedelta(days=i)
        date_str = str(day)
        date_list.append(date_str)
        date_list.append(date_str)
    return date_list