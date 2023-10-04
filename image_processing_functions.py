import numpy as np
from numpy.lib.stride_tricks import as_strided


def strided_rescale(g, bin_fac):
    strided = as_strided(g,
                         shape=(g.shape[0] // bin_fac, g.shape[1] // bin_fac, bin_fac, bin_fac),
                         strides=((g.strides[0] * bin_fac, g.strides[1] * bin_fac) + g.strides))
    return strided.mean(axis=-1).mean(axis=-1)


def convert_to_3D_array(array1):
    array1 = np.delete(array1, 2, axis=0)
    array1 = np.delete(array1, 2, axis=1)
    # array1 = np.delete(array1, 3, axis=2)

    return array1


def mask_saturation(rast, nodata=-9999):
    '''
    Masks out saturated values (surface reflectances of 16000). Arguments:
        rast    A gdal.Dataset or NumPy array
        nodata  The NoData value; defaults to -9999.
    '''
    # Can accept either a gdal.Dataset or numpy.array instance
    if not isinstance(rast, np.ndarray):
        rast = rast.ReadAsArray()

    # Create a baseline "nothing is saturated in any band" raster
    mask = np.empty((1, rast.shape[1], rast.shape[2]))
    mask.fill(False)

    # Update the mask for saturation in any band
    for i in range(rast.shape[0]):
        np.logical_or(mask,
                      np.in1d(rast[i, ...].ravel(), (16000,)).reshape(rast[i, ...].shape),
                      out=mask)

    # Repeat the NoData value across the bands
    np.place(rast, mask.repeat(rast.shape[0], axis=0), (nodata,))
    return rast