import numpy as np
from astropy.visualization.mpl_normalize import ImageNormalize
from image_processing_functions import convert_to_3D_array, strided_rescale
from astropy.visualization import AsinhStretch
import cv2
from datetime import date, timedelta
from sunpy.visualization.colormaps import color_tables as ct
import astropy.units as u


def creating_1600_image_list(list_1600, x, y, exposure_list):
    peak_x = x
    peak_y = 4096 - y
    bruh = 300
    sdoaia171 = ct.aia_color_table(1600 * u.angstrom)
    new_list = []
    full_1600_image_list = []
    intensity_list = []
    contour_list = []

    for i in range(0, len(list_1600)):
        data = np.clip(list_1600[i], 0, 1500)

        # data = ImageNormalize(stretch=AsinhStretch(
        #    0.01)).__call__(data)
        print("data")
        print(data.shape)
        cropped_low_resolution_image_data = data[peak_x - bruh:peak_x + bruh,
                                            peak_y - bruh:peak_y + bruh]
        intensity_value = np.divide(cropped_low_resolution_image_data.sum(), exposure_list[i])
        intensity_list.append(intensity_value)
        cropped_low_resolution_image_data = ImageNormalize(stretch=AsinhStretch(
            0.01)).__call__(cropped_low_resolution_image_data)

        cropped_low_resolution_image_data = (np.uint8(sdoaia171(cropped_low_resolution_image_data) * 255))
        low_resolution_image_data = strided_rescale(data, 8)
        low_resolution_image_data = ImageNormalize(stretch=AsinhStretch(
            0.01)).__call__(low_resolution_image_data)
        low_resolution_image_data = (np.uint8(sdoaia171(low_resolution_image_data) * 255))

        # cropped_low_resolution_image_data = (np.uint8(sdoaia171(cropped_low_resolution_image_data) * 255))

        new_list.append(cropped_low_resolution_image_data)
        full_1600_image_list.append(low_resolution_image_data)
    final_1600_list = []
    base_image = new_list[0]
    print("Zero frame")
    zero_frame = np.array(convert_to_3D_array(base_image.astype(np.int16)))
    for k in range(1, len(new_list)):
        img = new_list[k]
        one_frame = np.array(convert_to_3D_array(img).astype(np.int16))
        # im = np.subtract(zero_frame,one_frame)
        im = np.absolute(one_frame - zero_frame)
        imgray = im.astype(np.uint8)

        imgray = cv2.cvtColor(imgray, cv2.COLOR_BGR2GRAY)

        print(imgray)
        ret, thresh = cv2.threshold(imgray, 80, 255, cv2.THRESH_BINARY)
        contours, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)

        cnts = sorted(contours, key=cv2.contourArea, reverse=True)[:5]
        number_of_contours = 4
        try:
            for j in range(0, number_of_contours):
                cv2.drawContours(img, cnts[j], -1, (0, 255, 0), 1)
        except IndexError:
            print("less than 4 contours")
            pass

        final_1600_list.append(img)
        contour_list.append(cnts)
        print("final 1600 list")
        print(img)
        print(final_1600_list)
    return final_1600_list, full_1600_image_list, intensity_list, contour_list


def create_difference_and_ratio_lists(image_data_list, sdoaia):
    run_diff_list = []
    base_diff_list = []
    run_ratio_list = []
    base_ratio_list = []
    zero_image = image_data_list[0]
    for i in range(1, len(image_data_list)):
        first_image = image_data_list[i]
        second_image = image_data_list[i - 1]
        print(first_image)
        print(second_image)
        print("imagess")
        run_diff_image = strided_rescale(np.subtract(first_image, second_image), 8)
        run_diff_image = np.clip(run_diff_image, 0, 1500)
        run_diff_image = ImageNormalize(stretch=AsinhStretch(
            0.01)).__call__(run_diff_image)
        run_diff_image = (np.uint8(sdoaia(run_diff_image) * 255))

        base_diff_image = strided_rescale(np.subtract(first_image, zero_image), 8)
        base_diff_image = np.clip(base_diff_image, 0, 1500)
        base_diff_image = ImageNormalize(stretch=AsinhStretch(
            0.01)).__call__(base_diff_image)
        base_diff_image = (np.uint8(sdoaia(base_diff_image) * 255))

        run_ratio_image = strided_rescale(np.divide(first_image, second_image), 8)
        run_ratio_image = np.clip(run_ratio_image, 1, 10)
        run_ratio_image = ImageNormalize(stretch=AsinhStretch(
            0.01)).__call__(run_ratio_image)
        run_ratio_image = (np.uint8(sdoaia(run_ratio_image) * 255))

        base_ratio_image = strided_rescale(np.divide(first_image, zero_image), 8)
        base_ratio_image = np.clip(base_ratio_image, 1, 10)
        base_ratio_image = ImageNormalize(stretch=AsinhStretch(
            0.01)).__call__(base_ratio_image)
        base_ratio_image = (np.uint8(sdoaia(base_ratio_image) * 255))

        run_diff_list.append(run_diff_image)
        base_diff_list.append(base_diff_image)
        run_ratio_list.append(run_ratio_image)
        base_ratio_list.append(base_ratio_image)

    return run_diff_list, base_diff_list, run_ratio_list, base_ratio_list


def create_date_list(sdate, edate):
    print("fdsfdsfdsf")
    sdate = date(int(sdate[:4]), int(sdate[5:7]), int(sdate[8:10]))  # start date
    edate = date(int(edate[:4]), int(edate[5:7]), int(edate[8:10]))  # end date

    delta = edate - sdate
    date_list = []

    for i in range(delta.days):
        day = sdate + timedelta(days=i)

        date_list.append(str(day))
        date_list.append(str(day))
    return date_list