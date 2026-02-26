#!/usr/bin/env python3

# --- Aim
# The H&E images are very large and generally contain 2 tissues
# Crop the images to contain only the tissue of interest so they can be loaded more easily


# --- Define environment
# env: xen-hev

# --- Load packages
import argparse
import os
import subprocess

import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad

import tifffile


# --- Argparse arguments

parser = argparse.ArgumentParser(description="Import Visium data & save as .h5ad file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required
parser.add_argument("--raw_data_dir", help="Path to raw data")
parser.add_argument("--h5ad_dir", help="Path to directory for .h5ad output files")
parser.add_argument("--sample_id", help="Sample ID")
parser.add_argument("--project_dir", help="Path to the root of the project directory")
parser.add_argument("--outs_path", help="Path to spaceranger 'outs'")
parser.add_argument("--use_data", help="Which data to use (log_norm or SCT)")


args = vars(parser.parse_args())


# --- Arguments
raw_data_dir = args["raw_data_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]

outs_path=args["outs_path"]
project_dir=args['project_dir']


# --- Load data

# Load the file that has the pixel coordinates
adata = ad.read_h5ad(os.path.join(h5ad_dir, 'with_images', f'{sample_id}_with_image.h5ad'))


# --- Calculate crop coordinates based on the spots
# To automatically crop the image to where the data is

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

coords = auto_crop_xenium(adata, border_size=0.1, spatial_key='spatial')
print(f"Crop coordinates: {coords}")


# --- Crop the H&E image

image_dir = os.path.join(project_dir, 'data', 'images', f'{sample_id}')

# Find the .ome.tif file in the directory
image_file = [f for f in os.listdir(image_dir) if f.endswith('.ome.tif')]
if len(image_file) == 0:
    raise FileNotFoundError(f"No .ome.tif file found in {image_dir}")
elif len(image_file) > 1:
    raise ValueError(f"Multiple .ome.tif files found in {image_dir}: {image_file}")
else:
    image_file = image_file[0]

print(f"Using image file: {image_file}")
image_path = os.path.join(image_dir, image_file)

# Use tifffile to load the image
with tifffile.TiffFile(image_path) as tif:
    original_ome_metadata = tif.pages[0].description
    memmapped_array = tif.asarray()
    # Extract the crop coordinates
    x_start, x_end, y_start, y_end = coords
    # Ensure the coordinates are within the image bounds
    x_start = max(0, x_start)
    y_start = max(0, y_start)
    x_end = min(memmapped_array.shape[1], x_end)
    y_end = min(memmapped_array.shape[0], y_end)
    print(f"Cropping image to: x({x_start}:{x_end}), y({y_start}:{y_end})")
    print(f"Cropped image shape will be: ({y_end - y_start}, {x_end - x_start}, {memmapped_array.shape[2]})")
    # Perform the crop using memory mapping
    cropped_data = memmapped_array[y_start:y_end, x_start:x_end]

# Save the cropped data as a new OME-TIFF
filename = f'{sample_id}_cropped_HE'
os.makedirs(f'{image_dir}/cropped_images', exist_ok=True)

tifffile.imwrite(
    f'{image_dir}/cropped_images/{filename}.tif',
    cropped_data,
)

# Load John's QuPath module
# -- Optional: convert the cropped .tif to pyramidal .ome.tif format
# Instructions here: https://www.10xgenomics.com/support/software/xenium-explorer/latest/tutorials/xe-image-file-conversion

# cmd = (
#     "export MODULEPATH=/g/data/a56/software/modules:$MODULEPATH; "
#     "module load QuPath; "
#     f"cd {image_dir}/cropped_images; "
#     f"ndpi2ome {filename}.tif {filename}.ome.tif"
# )

# subprocess.run(cmd, shell=True, check=True)


# --- Now also crop the hires & lowres images in the object

hires_img = adata.uns['spatial'][sample_id]['images']['hires']
lowres_img = adata.uns['spatial'][sample_id]['images']['lowres']

scalefactor_hires = adata.uns['spatial'][sample_id]['scalefactors']['tissue_hires_scalef']
scalefactor_lowres = adata.uns['spatial'][sample_id]['scalefactors']['tissue_lowres_scalef']

# Calculate the crop coordinates for both images
x_start_hires = int(x_start * scalefactor_hires)
x_end_hires = int(x_end * scalefactor_hires)
y_start_hires = int(y_start * scalefactor_hires)
y_end_hires = int(y_end * scalefactor_hires)

x_start_lowres = int(x_start * scalefactor_lowres)
x_end_lowres = int(x_end * scalefactor_lowres)
y_start_lowres = int(y_start * scalefactor_lowres)
y_end_lowres = int(y_end * scalefactor_lowres)

# Crop the images
cropped_hires = hires_img[y_start_hires:y_end_hires, x_start_hires:x_end_hires]
cropped_lowres = lowres_img[y_start_lowres:y_end_lowres, x_start_lowres:x_end_lowres]


# --- Shift the coordinate system according to the crop coords

pixel_coords = adata.obsm['spatial'].copy()
pixel_coords[:, 0] = pixel_coords[:, 0] - x_start
pixel_coords[:, 1] = pixel_coords[:, 1] - y_start

# Update the obsm with the new coordinates
adata.obsm['spatial'] = pixel_coords

# Update the images in the object
adata.uns['spatial'][sample_id]['images']['hires'] = cropped_hires
adata.uns['spatial'][sample_id]['images']['lowres'] = cropped_lowres


# --- Save the adata object

adata.write_h5ad(os.path.join(h5ad_dir, 'with_images', f'{sample_id}_with_image.h5ad'))

sc.pl.spatial(adata, color='x_centroid', img_key='hires',)

# Also save the coords that were used for cropping
coords_df = pd.DataFrame([[x_start, x_end, y_start, y_end]], columns=['x_start', 'x_end', 'y_start', 'y_end'])
coords_df.to_csv(os.path.join(image_dir, 'cropped_images', f'{sample_id}_crop_coordinates.csv'), index=False)


# --- Finish

print('Script successfully completed!')