#!/usr/bin/env python3

# --- Aim
# Unpack the 10x pipeline outputs deposited on GEO (GSE351850) into the Space Ranger / Xenium Onboard
# Analysis folder layout that the other preprocessing scripts read:
#   data/processed/spaceranger_outs/visium/<section_id>/
#   data/processed/spaceranger_outs/visium_hd/<section_id>/
#   data/processed/xenium_outs/<section_id>/
# GEO stores each sample's files flat, gzipped and prefixed with the GSM accession and section ID
# (e.g. GSM10145003_P12_LN_Vis_spatial_tissue_positions.csv.gz). This script puts them back in their
# folders and decompresses them.
# GEO holds the Visium HD 2 um and 8 um bins only. The 16 um and 48 um bins are rebuilt from
# hd_feature_slice.h5 the way Space Ranger makes them: 2 um counts summed in 8x8 and 24x24 blocks, a bin
# is under tissue if any of its 2 um bins is, and bin centres are placed with the slide transform.

# --- Define environment
# env: vis-hd

# --- Load packages
import argparse
import glob
import gzip
import json
import os
import re
import shutil
import tarfile

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

# --- Argparse arguments

parser = argparse.ArgumentParser(description="Unpack the GEO deposit into the Space Ranger / Xenium Onboard Analysis folder layout",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--project_dir", required=True, help="Path to the root of the project directory")
parser.add_argument("--modality", required=True, choices=["visium", "visium_hd", "xenium"], help="Which platform to unpack")
parser.add_argument("--geo_dir", default=None,
                    help="Folder with the GEO files: the sub-series RAW.tar and/or loose GSM files (default: data/geo/<modality>)")
parser.add_argument("--out_dir", default=None, help="Root of the unpacked layout (default: data/processed)")
parser.add_argument("--samples", nargs="*", default=None, help="Section IDs to unpack (default: all in the sample sheet)")
parser.add_argument("--overwrite", action="store_true", help="Replace files that already exist")

args = vars(parser.parse_args())


# --- Arguments

project_dir = args["project_dir"]
modality = args["modality"]
geo_dir = args["geo_dir"] or os.path.join(project_dir, "data", "geo", modality)
out_root = args["out_dir"] or os.path.join(project_dir, "data", "processed")
overwrite = args["overwrite"]

out_dir = os.path.join(out_root, {"visium": "spaceranger_outs/visium",
                                  "visium_hd": "spaceranger_outs/visium_hd",
                                  "xenium": "xenium_outs"}[modality])


# --- Section IDs from the sample sheets

sheets = {"visium": ["sample_sheet_visium.csv", "sample_sheet_visium_lymph_nodes.csv"],
          "visium_hd": ["sample_sheet_visium_hd.csv"],
          "xenium": ["sample_sheet_xenium.csv"]}[modality]
section_ids = pd.concat([pd.read_csv(os.path.join(project_dir, "config", s)) for s in sheets])["section_id"].tolist()
if args["samples"]:
    section_ids = [s for s in section_ids if s in args["samples"]]
# Longest first, so P12_LN_Vis is matched before a shorter ID that it starts with
section_ids = sorted(section_ids, key=len, reverse=True)
print(f"Unpacking {len(section_ids)} {modality} sections from {geo_dir} into {out_dir}")


# --- Map a GEO file name to its place in the 10x layout

def target_path(geo_name):
    """GSM<digits>_<section_id>_<name> -> (section_id, path inside the section folder), or None."""
    m = re.match(r"GSM\d+_(.+)$", geo_name)
    if not m:
        return None
    rest = m.group(1)
    section = next((s for s in section_ids if rest.startswith(s + "_")), None)
    if section is None:
        return None
    name = rest[len(section) + 1:]
    if name.endswith(".gz"):
        name = name[:-3]                       # GEO's compression, not part of the 10x file name
    if modality == "xenium":
        if name.startswith("morphology_focus_morphology_focus_"):
            name = "morphology_focus/" + name[len("morphology_focus_"):]
        return section, name
    prefix = ""
    binned = re.match(r"(square_\d{3}um)_(.+)$", name)
    if binned:
        prefix, name = f"binned_outputs/{binned.group(1)}/", binned.group(2)
    if name.startswith("spatial_"):
        name = "spatial/" + name[len("spatial_"):]
    return section, prefix + name


def write_member(fileobj, geo_name, dest):
    """Stream one GEO file to dest, decompressing it if GEO gzipped it."""
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    src = gzip.GzipFile(fileobj=fileobj) if geo_name.endswith(".gz") else fileobj
    with open(tmp, "wb") as out:
        shutil.copyfileobj(src, out, length=16 * 1024 * 1024)
    os.replace(tmp, dest)


# --- Unpack the RAW.tar members and any loose GSM files

written = skipped = 0


def handle(geo_name, open_member):
    global written, skipped
    target = target_path(os.path.basename(geo_name))
    if target is None:
        return
    section, rel = target
    dest = os.path.join(out_dir, section, rel)
    if os.path.exists(dest) and not overwrite:
        skipped += 1
        return
    with open_member() as f:
        write_member(f, geo_name, dest)
    written += 1


for tar_path in sorted(glob.glob(os.path.join(geo_dir, "*_RAW.tar"))):
    print(f"Reading {os.path.basename(tar_path)}")
    with tarfile.open(tar_path, "r:") as tar:
        for member in tar:
            if member.isfile():
                handle(member.name, lambda m=member: tar.extractfile(m))

for loose in sorted(glob.glob(os.path.join(geo_dir, "**", "GSM*"), recursive=True)):
    if os.path.isfile(loose):
        handle(loose, lambda p=loose: open(p, "rb"))

print(f"Files written: {written}; already present (skipped): {skipped}")


# --- Visium HD: image links per bin size and the 16 um / 48 um bins

def read_10x_h5(path):
    """10x HDF5 matrix -> (CSC genes x barcodes, barcodes)."""
    with h5py.File(path, "r") as f:
        m = f["matrix"]
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(m["shape"][:]))
        barcodes = m["barcodes"][:]
    return X, barcodes


def write_10x_h5(path, X, barcodes, template):
    """Write a 10x HDF5 matrix (genes x barcodes, CSC), copying features and file attributes from template."""
    tmp = path + ".part"
    with h5py.File(template, "r") as src, h5py.File(tmp, "w") as dst:
        for key, value in src.attrs.items():
            dst.attrs[key] = value
        m = dst.create_group("matrix")
        src.copy(src["matrix/features"], m, name="features")
        m.create_dataset("barcodes", data=np.asarray(barcodes, dtype=src["matrix/barcodes"].dtype))
        m.create_dataset("data", data=X.data.astype(src["matrix/data"].dtype), compression="gzip")
        m.create_dataset("indices", data=X.indices.astype(src["matrix/indices"].dtype), compression="gzip")
        m.create_dataset("indptr", data=X.indptr.astype(src["matrix/indptr"].dtype), compression="gzip")
        m.create_dataset("shape", data=np.array(X.shape, dtype=src["matrix/shape"].dtype))
    os.replace(tmp, path)


def rebin(section_dir, k):
    """Build binned_outputs/square_<2k>um from hd_feature_slice.h5 (k = 2 um bins per side)."""
    size = f"square_{2 * k:03d}um"
    bin_dir = os.path.join(section_dir, "binned_outputs", size)
    matrix_file = os.path.join(bin_dir, "filtered_feature_bc_matrix.h5")
    if os.path.exists(matrix_file) and not overwrite:
        return
    os.makedirs(os.path.join(bin_dir, "spatial"), exist_ok=True)
    template = os.path.join(section_dir, "binned_outputs", "square_008um", "filtered_feature_bc_matrix.h5")
    with h5py.File(template, "r") as f:
        matrix_ids = f["matrix/features/id"][:]
    n_genes = len(matrix_ids)

    with h5py.File(os.path.join(section_dir, "hd_feature_slice.h5"), "r") as fs:
        meta = json.loads(fs.attrs["metadata_json"])
        n_rows, n_cols = int(meta["nrows"]), int(meta["ncols"])
        transform = np.array(meta["transform_matrices"]["spot_colrow_to_microscope_colrow"]).reshape(3, 3)
        grid_rows, grid_cols = -(-n_rows // k), -(-n_cols // k)
        # The feature slices are numbered by the feature slice's own feature list; map them to the matrix rows
        row_of = pd.Series(np.arange(n_genes), index=matrix_ids)
        slice_ids = fs["features/id"][:]
        # Raw 2 um counts of every gene, all spots (including those outside the tissue), summed per block
        genes, blocks, counts = [], [], []
        for key, g in fs["feature_slices"].items():
            feature_id = slice_ids[int(key)]
            if feature_id not in row_of.index:
                continue
            r = g["row"][:].astype(np.int64)
            c = g["col"][:].astype(np.int64)
            genes.append(np.full(r.size, row_of[feature_id], dtype=np.int64))
            blocks.append((r // k) * grid_cols + (c // k))
            counts.append(g["data"][:].astype(np.int64))
    # Bins under tissue: those with any 2 um bin under tissue
    spots = pd.read_parquet(os.path.join(section_dir, "binned_outputs", "square_002um", "spatial", "tissue_positions.parquet"),
                            columns=["in_tissue", "array_row", "array_col"])
    spots = spots[spots["in_tissue"] == 1]
    in_tissue = np.unique((spots["array_row"].to_numpy(np.int64) // k) * grid_cols + spots["array_col"].to_numpy(np.int64) // k)
    X = sp.csc_matrix((np.concatenate(counts), (np.concatenate(genes), np.concatenate(blocks))),
                      shape=(n_genes, grid_rows * grid_cols))
    X = X[:, in_tissue]
    X.sort_indices()

    def barcode(idx):
        return [f"s_{2 * k:03d}um_{i // grid_cols:05d}_{i % grid_cols:05d}-1" for i in idx]

    write_10x_h5(matrix_file, X, barcode(in_tissue), template)

    # Positions of all bins: the centre of each block, mapped through the slide transform
    all_idx = np.arange(grid_rows * grid_cols)
    array_row, array_col = all_idx // grid_cols, all_idx % grid_cols
    centre = (k - 1) / 2
    xyw = transform @ np.vstack([array_col * k + centre, array_row * k + centre, np.ones(all_idx.size)])
    reference = pd.read_parquet(os.path.join(section_dir, "binned_outputs", "square_008um", "spatial", "tissue_positions.parquet"))
    positions = pd.DataFrame({
        "barcode": barcode(all_idx),
        "in_tissue": np.isin(all_idx, in_tissue).astype(reference["in_tissue"].dtype),
        "array_row": array_row.astype(reference["array_row"].dtype),
        "array_col": array_col.astype(reference["array_col"].dtype),
        "pxl_row_in_fullres": xyw[1] / xyw[2],
        "pxl_col_in_fullres": xyw[0] / xyw[2],
    })
    positions.to_parquet(os.path.join(bin_dir, "spatial", "tissue_positions.parquet"), index=False)

    # Scale factors: the 2 um ones, with the bin size and spot diameter scaled up
    with open(os.path.join(section_dir, "binned_outputs", "square_002um", "spatial", "scalefactors_json.json")) as f:
        scalefactors = json.load(f)
    scalefactors["spot_diameter_fullres"] = scalefactors["spot_diameter_fullres"] * k
    scalefactors["bin_size_um"] = float(2 * k)
    with open(os.path.join(bin_dir, "spatial", "scalefactors_json.json"), "w") as f:
        json.dump(scalefactors, f, indent=4)
    print(f"  {os.path.basename(section_dir)}: rebuilt {size} ({X.shape[1]} bins under tissue)")


if modality == "visium_hd":
    for section in section_ids:
        section_dir = os.path.join(out_dir, section)
        if not os.path.isdir(section_dir):
            continue
        for k in (8, 24):
            rebin(section_dir, k)
        # Space Ranger keeps a copy of the slide images in each bin size's spatial/ folder
        for bin_dir in glob.glob(os.path.join(section_dir, "binned_outputs", "square_*um")):
            for image in glob.glob(os.path.join(section_dir, "spatial", "*")):
                link = os.path.join(bin_dir, "spatial", os.path.basename(image))
                if not os.path.lexists(link):
                    os.symlink(os.path.relpath(image, os.path.dirname(link)), link)

print("Script successfully completed")
