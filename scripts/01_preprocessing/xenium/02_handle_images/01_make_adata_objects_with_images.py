#!/usr/bin/env python3

# --- Aim
# Use alignment matrix to transform the Xenium coordinates to aligned pixel coordinates in image space
# Then add the image to the AnnData object and save it


# --- Define environment
# env: tifffile

# --- Load packages
import argparse
import os

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
parser.add_argument("--min_count", help="Minimal number of counts per spot", type=int)
parser.add_argument("--min_gene", help="Minimal number of genes per spot", type=int)
parser.add_argument("--outs_path", help="Path to spaceranger 'outs'")
parser.add_argument("--use_data", help="Which data to use (log_norm or SCT)")


args = vars(parser.parse_args())


# --- Arguments
raw_data_dir = args["raw_data_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]

min_count=args["min_count"]
min_gene=args["min_gene"]

outs_path=args["outs_path"]
project_dir=args['project_dir']


# --- Load the sample sheet

sample_sheet = pd.read_csv(os.path.join(project_dir, 'config', 'sample_sheet.csv'))

# Extract the slide id from the sample sheet
xen_slide_id = sample_sheet.loc[sample_sheet['section_id'] == sample_id, 'xen_slide_id'].values[0]
print(f"Xenium slide ID: {xen_slide_id}")


# --- Load data

# Load the log_norm file
adata = ad.read_h5ad(os.path.join(h5ad_dir, 'log_norm', f'{sample_id}_log_norm.h5ad'))


# --- Load the transformation matrix

alignment_dir = os.path.join(project_dir, 'resources', 'image_alignment', '5k_alignment', f'{sample_id}_alignment_files')
transform_matrix = pd.read_csv(os.path.join(alignment_dir, 'matrix.csv'), header=None).to_numpy()


# --- Try to transform the coordinates of anndata by doing an inverse transform

# Calculate the inverse transformation matrix
inverse_transform = np.linalg.inv(transform_matrix)

# Extract the Xenium coordinates from the AnnData object
xenium_coords = adata.obsm['spatial']

# First convert from microns to pixels
microns_per_pixel = 0.2125  # Assuming 0.2125 microns per pixel, can be found in experiment.xenium file
xenium_coords = xenium_coords / microns_per_pixel  # Convert to pixels

# Prepare coordinates for matrix multiplication
#    We need to add a "1" to each coordinate (homogeneous coordinate)
#    to work with the 3x3 transformation matrix.
#    This creates a new array of shape (n_cells, 3)
n_cells = xenium_coords.shape[0]
hom_coords = np.hstack([xenium_coords, np.ones((n_cells, 1))])

# Apply the inverse transformation
transformed_coords = (inverse_transform @ hom_coords.T).T

# Extract the pixel coords
pixel_coords = transformed_coords[:, :2]

# Update the AnnData object with the new pixel coordinates
adata.obsm['spatial'] = pixel_coords

# Save xenium coords as well
adata.obsm['xenium_um_spatial'] = xenium_coords


# --- Load the image

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

# Inspect the pyramid levels in the OME-TIFF
with tifffile.TiffFile(image_path) as tif:
    # A TIFF file can contain multiple data series (e.g., time series).
    # For a standard pyramidal image, we are interested in the first series.
    pyramid = tif.series[0]
    print(f"Inspecting pyramid in: '{image_path}'")
    print(f"Found {len(pyramid.levels)} levels in the pyramid.")
    print("-" * 30)
    # Iterate through each level in the pyramid series
    for i, level in enumerate(pyramid.levels):
        # The .shape attribute contains the dimensions (height, width, channels, etc.)
        # This reads metadata only, not the pixel data itself.
        print(f"Level {i}: Dimensions = {level.shape}, Data Type = {level.dtype}")

# Load level 4 & the lowest level as hires & lowres
# Load lowest resolution level of the pyramid
with tifffile.TiffFile(image_path) as tif:
    pyramid = tif.series[0]
    lowest_res_level = pyramid.levels[-1]  # Last level is the lowest resolution
    lowres_img = lowest_res_level.asarray()  # This reads the pixel data into memory
    hires_level = pyramid.levels[4]  # First level is the highest resolution
    hires_img = hires_level.asarray()  # This reads the pixel data into memory

# Calculate a scalefactor between for both
scalefactor_lowres = pyramid.levels[0].shape[0] / lowest_res_level.shape[0]
scalefactor_hires = pyramid.levels[0].shape[0] / hires_level.shape[0]

# Calculate a sensible spot_size
spot_size = 40/microns_per_pixel

# Add img & scalefactor to the AnnData object
adata.uns['spatial'] = {sample_id: {
        'images': {
            'hires': hires_img,
            'lowres': lowres_img
        },
        'scalefactors': {
            'tissue_hires_scalef': 1.0 / (scalefactor_hires),
            'tissue_lowres_scalef': 1.0 / (scalefactor_lowres),
            'spot_diameter_fullres': spot_size,
            'microns_per_pixel': microns_per_pixel,
        }
    }
}


# --- Save the adata object with the image & transformed coordinates

out_dir = os.path.join(h5ad_dir, 'with_images')
os.makedirs(out_dir, exist_ok=True)

adata.write_h5ad(os.path.join(out_dir, f'{sample_id}_with_image.h5ad'))


# --- Finish

print('Script successfully completed')
