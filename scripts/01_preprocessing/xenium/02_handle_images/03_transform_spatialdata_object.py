#!/usr/bin/env python3

# --- Aim
# Transform the spatialdata objects so they are aligned with the H&E image
# This will make it easier to plot zoomed in images with the segmentation mask and transcripts on top


# --- Define environment
# env: spatdat

# --- Load packages
import argparse
import os
import shutil

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import spatialdata as sd
from spatialdata.transformations import Affine, Translation, Sequence, Scale, set_transformation


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
project_dir = args["project_dir"]
raw_data_dir = args["raw_data_dir"]
h5ad_dir=args["h5ad_dir"]
sample_id=args["sample_id"]

min_count=args["min_count"]
min_gene=args["min_gene"]

outs_path=args["outs_path"]

# Make directories if they don't exist yet
os.makedirs(h5ad_dir, exist_ok=True)


# --- Load the sample sheet

sample_sheet = pd.read_csv(os.path.join(project_dir, 'config', 'sample_sheet.csv'))

# Extract the slide id from the sample sheet
xen_slide_id = sample_sheet.loc[sample_sheet['section_id'] == sample_id, 'xen_slide_id'].values[0]
print(f"Xenium slide ID: {xen_slide_id}")


# --- Load the previously created spatialdata object

spat_dat_dir = os.path.join(project_dir, 'data', 'processed', 'spatialdata_objects', 'log_norm')

filename = f'{sample_id}.zarr'

# Might need to unzip first
if not os.path.exists(os.path.join(spat_dat_dir, filename)):
    print(f"Unzipping {filename}")
    shutil.unpack_archive(os.path.join(spat_dat_dir, filename.replace('zarr', 'zip')), os.path.join(spat_dat_dir, filename))

sdata = sd.read_zarr(os.path.join(spat_dat_dir, filename))


# --- Load the transformation matrix
# Records how the image was aligned to the xenium coordinates

alignment_dir = os.path.join(project_dir, 'resources', 'image_alignment', '5k_alignment', f'{sample_id}_alignment_files')
transform_matrix = pd.read_csv(os.path.join(alignment_dir, 'matrix.csv'), header=None).to_numpy()


# --- Do an inverse transform to align the data with the image

# Calculate the inverse transformation matrix
inv_transform_matrix = np.linalg.inv(transform_matrix)

# Define the affine transformation
affine_transformation = Affine(inv_transform_matrix, input_axes=('x', 'y'), output_axes=('x', 'y'))


# --- Because the image was cropped, we need to do an additional shift transformation

# Load the crop coords
image_dir = os.path.join(project_dir, 'data', 'images', sample_id, 'cropped_images')
coords_file = f'{sample_id}_crop_coordinates.csv'

coords_df = pd.read_csv(os.path.join(image_dir, coords_file))
coords = coords_df.iloc[0].to_list()
x_start, x_end, y_start, y_end = coords

# Create the separate shift transformation
shift_transformation = Translation([-y_start, -x_start], axes=('y', 'x'))

# Microns per pixel
microns_per_pixel = 0.2125  # Assuming 0.2125 microns per pixel, can be found in experiment.xenium file

# The xenium coordinates are in microns, the image in pixels, so we need to scale as well
# Define scaling to go from microns to pixels
scaling = Scale([1/microns_per_pixel, 1/microns_per_pixel], axes=('x', 'y'))

# Combine the affine and shift transformations into a sequence
combined_transformation = Sequence([scaling, affine_transformation, shift_transformation])


# --- Apply the transformation to all elements in the SpatialData object

# Apply the single combined transformation to all elements
for key, shapes in sdata.shapes.items():
    set_transformation(shapes, combined_transformation, to_coordinate_system='image_space')

for key, points in sdata.points.items():
    set_transformation(points, combined_transformation, to_coordinate_system='image_space')

for key, labels in sdata.labels.items():
    set_transformation(labels, combined_transformation, to_coordinate_system='image_space')


# The object now contains the transformation info, but the coordinates are not yet transformed
# To actually transform the coordinates, use sd.transform()
# This will create a new spatialdata object with the transformed coordinates
# Better to do this on the fly because it throws an error when trying to save the object otherwise
# transformed_sdata = sd.transform(sdata, to_coordinate_system='image_space')


# --- Save the spatialdata object with new alignment

# save AnnData object to *.h5ad file
outs_dir = os.path.join(project_dir, 'data', 'processed', 'spatialdata_objects', 'aligned_with_HE')
os.makedirs(outs_dir, exist_ok=True)

print(f"Saving .h5ad files to {outs_dir}")

# Save the sdata object that contains the transformation info
sdata.write(os.path.join(outs_dir, f"{sample_id}.zarr"))

# Zip the folder to avoid having thousands of file in the dvc object
shutil.make_archive(os.path.join(outs_dir, f"{sample_id}"), 'zip', os.path.join(outs_dir, f"{sample_id}.zarr"))

# Remove the unzipped folder to save space
shutil.rmtree(os.path.join(outs_dir, f"{sample_id}.zarr"))

# Also remove the unzipped folder from the previous step if it still exists, as well as the previous object
# Just because the objects are quite large and this one contains the same info plus the transformation info
if os.path.exists(os.path.join(spat_dat_dir, filename)):
    shutil.rmtree(os.path.join(spat_dat_dir, filename))

if os.path.exists(os.path.join(spat_dat_dir, filename.replace('zarr', 'zip'))):
    os.remove(os.path.join(spat_dat_dir, filename.replace('zarr', 'zip')))

print("Script succesfully completed")