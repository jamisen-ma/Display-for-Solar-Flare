"""Image processing helpers."""

from __future__ import annotations

import numpy as np
from numpy.lib.stride_tricks import as_strided
from typing import Union


def strided_rescale(image: np.ndarray, bin_factor: int) -> np.ndarray:
    """Downsample a 2D array by averaging non-overlapping ``bin_factor`` blocks."""
    if bin_factor <= 0:
        raise ValueError("bin_factor must be positive")
    h, w = image.shape[:2]
    if h < bin_factor or w < bin_factor:
        raise ValueError("image smaller than bin_factor")
    new_h = h // bin_factor
    new_w = w // bin_factor
    strided = as_strided(
        image,
        shape=(new_h, new_w, bin_factor, bin_factor),
        strides=((image.strides[0] * bin_factor, image.strides[1] * bin_factor) + image.strides[:2]),
    )
    return strided.mean(axis=-1).mean(axis=-1)


def convert_to_3D_array(array1: np.ndarray) -> np.ndarray:
    """Trim a 2D or 3D array to 3D with at least 2x2 spatial dims by removing last row/col.

    This mirrors the original behavior which deleted index 2 along first two axes.
    """
    if array1.shape[0] > 2:
        array1 = np.delete(array1, 2, axis=0)
    if array1.shape[1] > 2:
        array1 = np.delete(array1, 2, axis=1)
    return array1


def mask_saturation(rast: Union[np.ndarray, "gdal.Dataset"], nodata: float = -9999) -> np.ndarray:  # type: ignore[name-defined]
    """Mask saturated values (exactly 16000) with ``nodata`` across all bands."""
    if not isinstance(rast, np.ndarray):
        rast = rast.ReadAsArray()

    mask = np.zeros((1, rast.shape[1], rast.shape[2]), dtype=bool)
    for i in range(rast.shape[0]):
        saturated = np.in1d(rast[i, ...].ravel(), (16000,)).reshape(rast[i, ...].shape)
        np.logical_or(mask, saturated, out=mask)

    np.place(rast, mask.repeat(rast.shape[0], axis=0), nodata)
    return rast