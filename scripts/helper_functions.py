# Imports

import colorsys
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# --- Automatic crop

def auto_crop(adata, border_size=0.2):
    """
    Calculate crop coordinates to pass to squidpy.pl.spatial_scatter().
    
    :adata: An AnnData object containing Visium data
    
    :border_size: Size of the border around the spots
    """

    # get the coordinates of the spots
    coord = adata.obsm['spatial']

    coord_x = coord[:, 0]
    coord_y = coord[:, 1]

    # do some calculations
    x_min = np.min(coord_x)
    x_max = np.max(coord_x)

    y_min = np.min(coord_y)
    y_max = np.max(coord_y)

    x_range = x_max - x_min
    y_range = y_max - y_min

    # define a mid point
    mid_point = np.array([x_range/2 + x_min, y_range/2 + y_min])
    mid_point = mid_point.astype(int)

    if x_range > y_range:
        mid_dist = int(x_range/2)
        border = int((border_size/2)*x_range)
    else:
        mid_dist = int(y_range/2)
        border = int((border_size/2)*y_range)

    # define the crop coordinates for the square
    coord_square = [mid_point[0] - mid_dist - border, mid_point[1] - mid_dist -
                    border, mid_point[0] + mid_dist + border, mid_point[1] + mid_dist + border]

    # need to check that the coordinates are not out of bounds
    sample_id = list(adata.uns['spatial'].keys())[0]
    hires = adata.uns['spatial'][sample_id]['images']['hires']

    scale_factor = adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']

    img_max_y = list(hires.shape)[0] / scale_factor # find the max image dimensions
    img_max_x = list(hires.shape)[1] / scale_factor

    img_max_y = int(img_max_y)
    img_max_x = int(img_max_x)

    # shift the square if it goes out of bounds
    if coord_square[0] < 0:
        coord_square[2] -= coord_square[0]
        coord_square[0] = 0

    if coord_square[1] < 0:
        coord_square[3] -= coord_square[1]
        coord_square[1] = 0

    if coord_square[2] > img_max_x:
        coord_square[0] -= (img_max_x - coord_square[2])
        coord_square[2] = img_max_x

    if coord_square[3] > img_max_y:
        coord_square[1] -= (img_max_y - coord_square[3])
        coord_square[3] = img_max_y

    coord_square=tuple(coord_square)

    return coord_square

# The function above is less reliable when the hires image is smaller than the area where spots were captured
# Squidpy doesn't support cropping outside of the hires image area so that way you miss a lot of spots
# Therefore create an identical function that will work with scanpy, which does support cropping outside the image area

def auto_crop_scanpy(adata, border_size=0.05, spot_diameter=350):
    """
    Calculate crop coordinates to pass to scanpy.pl.scatter().
    
    :adata: An AnnData object containing Visium data
    
    :border_size: Size of the border around the spots

    :spot_diameter: Diameter of a Visium spot (only necessary if anndata contains just a single spot)
    """
    # get the coordinates of the spots
    coord = adata.obsm['spatial']
    coord_x = coord[:, 0]
    coord_y = coord[:, 1]
    # do some calculations
    x_min = np.min(coord_x)
    x_max = np.max(coord_x)
    y_min = np.min(coord_y)
    y_max = np.max(coord_y)
    x_range = x_max - x_min
    y_range = y_max - y_min
    # define a mid point
    mid_point = np.array([x_range / 2 + x_min, y_range / 2 + y_min])
    mid_point = mid_point.astype(int)
    if x_range == 0 and y_range == 0:
        # Only one spot, define mid_dist and border based on a fixed size
        mid_dist = spot_diameter
        border = int((border_size / 2) * mid_dist)
    else:
        if x_range > y_range:
            mid_dist = int(x_range / 2)
            border = int((border_size / 2) * x_range)
        else:
            mid_dist = int(y_range / 2)
            border = int((border_size / 2) * y_range)
    # define the crop coordinates for the square
    coord_square = [
        mid_point[0] - mid_dist - border,
        mid_point[0] + mid_dist + border,
        mid_point[1] - mid_dist - border,
        mid_point[1] + mid_dist + border,
    ]
    coord_square = tuple(coord_square)
    return coord_square


# Define a cropping function for Xenium - to make the plotting area square in the spatial plots
# Assuming coordinates are stored in adata.obsm['spatial]
# Output is a dictionary of coordinates that can be passed to crop_coord

def auto_crop_xenium(adata, border_size=0.1, spatial_key='spatial'):
    """
    Calculate crop coordinates to pass to sc.pl.spatial_scatter() for Xenium samples (make sure the plotting area is a square)
    Returns a dictionary of crop coordinates
    Pass the coordinates to the crop_coord parameter of scanpy's sc.pl.spatial_scatter()
    
    :adata: an AnnData object
    :border_size: size of the border around the cells (data points)
    :spatial_key: key in adata.obsm where the spatial coordinates are stored
    """ 

    # get the coordinates of the cells
    coord = adata.obsm[spatial_key]
    coord_x = coord[:, 0]
    coord_y = coord[:, 1]
    # do some calculations
    x_min = np.min(coord_x)
    x_max = np.max(coord_x)
    y_min = np.min(coord_y)
    y_max = np.max(coord_y)
    x_range = x_max - x_min
    y_range = y_max - y_min
    # define a mid point
    mid_point = np.array([x_range/2 + x_min, y_range/2 + y_min])
    mid_point = mid_point.astype(int)
    if x_range > y_range:
        mid_dist = int(x_range/2)
        border = int((border_size/2)*x_range)
    else:
        mid_dist = int(y_range/2)
        border = int((border_size/2)*y_range)
    # define the crop coordinates for the square
    coord_square = [mid_point[0] - mid_dist - border, mid_point[0] + mid_dist + border, mid_point[1] - mid_dist -
                border, mid_point[1] + mid_dist + border]
    coord_square = tuple(coord_square)

    return coord_square


# --- Color functions
# Description: functions for generating / altering colors

def generate_color_variations(hex_color, n, step_size=0.2):

    """
    Generate variations of a given color. Returns a list of darker and lighter variations of the given color, including the color itself.

    :hex_color: A color in hex format (e.g. #1f77b4)
    :n: The number of color variations to generate
    :step_size: Size factor for how much lighter / darker each shade should get (set higher for more distinct shades)
    """

    # Remove the #
    hex_color = hex_color.lstrip('#')
    # Convert hex color to RGB tuple
    r, g, b = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
    # Convert RGB values to HLS values
    h, l, s = colorsys.rgb_to_hls(r/255, g/255, b/255)
    # Generate n color variations
    colors = []

    if l > 0.5:
        # For light colors: 2/3 darker shades, 1/3 lighter shades
        # -1 because we're adding the original color back in the middle
        n_dark = math.ceil(n/3*2) - 1
        n_light = n - n_dark - 1
        for i in range(n_dark):
            # Adjust the lightness value to create lighter or darker shades
            l_new = max(min(l - (i+1)*step_size, 0.95), 0.05)
            # Convert the HLS values back to RGB values
            r_new, g_new, b_new = tuple(round(val*255) for val in colorsys.hls_to_rgb(h, l_new, s))
            # Convert the RGB values back to hex color
            hex_new = "#{:02x}{:02x}{:02x}".format(r_new, g_new, b_new)
            colors.append(hex_new)
        # Reverse the list so the darkest shade comes first
        colors.reverse()
        # Add the original color back in
        colors.append(f'#{hex_color}')
        # Then get some lighter shades
        for i in range(n_light):
            # Adjust the lightness value to create lighter or darker shades
            l_new = max(min(l + (i+1)*step_size, 0.95), 0.05)
            # Convert the HLS values back to RGB values
            r_new, g_new, b_new = tuple(round(val*255) for val in colorsys.hls_to_rgb(h, l_new, s))
            # Convert the RGB values back to hex color
            hex_new = "#{:02x}{:02x}{:02x}".format(r_new, g_new, b_new)
            colors.append(hex_new)

    else:
        # For darker colors: 1/3 darker shades, 2/3 lighter shades
        # -1 because we're adding the original color back in the middle
        n_dark = math.ceil(n/3*1) - 1
        n_light = n - n_dark - 1
        for i in range(n_dark):
            # Adjust the lightness value to create lighter or darker shades
            l_new = max(min(l - (i+1)*step_size, 0.95), 0.05)
            # Convert the HLS values back to RGB values
            r_new, g_new, b_new = tuple(round(val*255) for val in colorsys.hls_to_rgb(h, l_new, s))
            # Convert the RGB values back to hex color
            hex_new = "#{:02x}{:02x}{:02x}".format(r_new, g_new, b_new)
            colors.append(hex_new)
        # Reverse the list so the darkest shade comes first
        colors.reverse()
        # Add the original color back in
        colors.append(f'#{hex_color}')
        # Then get some lighter shades
        for i in range(n_light):
            # Adjust the lightness value to create lighter or darker shades
            l_new = max(min(l + (i+1)*step_size, 0.95), 0.05)
            # Convert the HLS values back to RGB values
            r_new, g_new, b_new = tuple(round(val*255) for val in colorsys.hls_to_rgb(h, l_new, s))
            # Convert the RGB values back to hex color
            hex_new = "#{:02x}{:02x}{:02x}".format(r_new, g_new, b_new)
            colors.append(hex_new)

    return colors


# --- Save plots as one pdf
# Note: not my code pulled this off google
# Plots all the plots of your script in a single pdf, each plot is a separate page in the pdf

def save_multi_image(filename):
    """
    Save all plots in a python script as a single pdf files.
    Each plot will have its own page.
    
    :filename: Name of the pdf to save the plots to ('example.pdf')
    :type filename: str
    :dpi: Resolution in dots per inch
    """

    pp = PdfPages(filename)
    fig_nums = plt.get_fignums()
    
    figs = [plt.figure(n) for n in fig_nums]
    for fig in figs:
        fig.savefig(pp, format='pdf', bbox_inches='tight')
    pp.close()


# --- Region reads from the aligned H&E OME-TIFFs
# The aligned images are whole slides (~10-13 GB in RAM at full resolution), so the zoom scripts
# read only the 1024 px tiles that overlap a crop, and overviews from a low-resolution pyramid
# level. Tiled reads via tifffile.aszarr need zarr>=3, which the spatdat env does not have, so
# DEFLATE tiles are decoded directly (anything else falls back to a full read).

def tiff_shape(path):
    """(height, width) of the full-resolution level."""
    import tifffile
    with tifffile.TiffFile(path) as tf:
        page = tf.pages[0]
        return int(page.imagelength), int(page.imagewidth)


def read_tiff_region(path, x0, x1, y0, y1):
    """Full-resolution RGB crop [y0:y1, x0:x1], clipped to the image like numpy slicing would be.
    Returns the crop and the clipped (x0, x1, y0, y1)."""
    import zlib
    import tifffile
    with tifffile.TiffFile(path) as tf:
        page = tf.pages[0]
        width, height, n_ch = int(page.imagewidth), int(page.imagelength), int(page.samplesperpixel)
        x0, x1 = max(int(x0), 0), min(int(x1), width)
        y0, y1 = max(int(y0), 0), min(int(y1), height)
        if x1 <= x0 or y1 <= y0:
            return np.zeros((max(y1 - y0, 0), max(x1 - x0, 0), 3), np.uint8), (x0, x1, y0, y1)
        fast =(page.is_tiled and int(page.predictor) == 1 and int(page.planarconfig or 1) == 1
                and int(page.compression) in (8, 32946) and page.dtype == np.uint8)
        if not fast:
            return page.asarray()[y0:y1, x0:x1, :3], (x0, x1, y0, y1)
        tw, tl = int(page.tilewidth), int(page.tilelength)
        tiles_across = (width + tw - 1) // tw
        tx0, tx1 = x0 // tw, (x1 - 1) // tw
        ty0, ty1 = y0 // tl, (y1 - 1) // tl
        canvas = np.zeros(((ty1 - ty0 + 1) * tl, (tx1 - tx0 + 1) * tw, n_ch), np.uint8)
        fh, offsets, counts = tf.filehandle, page.dataoffsets, page.databytecounts
        for ty in range(ty0, ty1 + 1):
            for tx in range(tx0, tx1 + 1):
                idx = ty * tiles_across + tx
                if counts[idx] == 0:
                    continue
                fh.seek(offsets[idx])
                tile = np.frombuffer(zlib.decompress(fh.read(counts[idx])), np.uint8).reshape(tl, tw, n_ch)
                canvas[(ty - ty0) * tl:(ty - ty0 + 1) * tl, (tx - tx0) * tw:(tx - tx0 + 1) * tw] = tile
        crop = canvas[y0 - ty0 * tl:y1 - ty0 * tl, x0 - tx0 * tw:x1 - tx0 * tw, :3]
    return crop, (x0, x1, y0, y1)


def read_tiff_overview(path, min_side=4000):
    """The coarsest pyramid level whose longer side is >= min_side px (RGB). Draw it with
    imshow(..., extent=(0, width, height, 0)) to keep full-resolution pixel coordinates."""
    import tifffile
    with tifffile.TiffFile(path) as tf:
        levels = tf.series[0].levels
        level = next((lv for lv in reversed(levels) if max(lv.shape[:2]) >= min_side), levels[0])
        return level.asarray()[..., :3]