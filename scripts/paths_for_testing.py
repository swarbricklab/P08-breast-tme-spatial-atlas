# This scripts contains paths for quick interactive testing of scripts on specific samples

import os
import pandas as pd


# --- Arguments
project_dir = os.getcwd()

data_dir = os.path.join(project_dir, 'data')
res_dir = os.path.join(project_dir, "results")
h5ad_dir = os.path.join(project_dir, "data/anndata_objects")
sample_sheet_path = os.path.join(project_dir, 'config', 'sample_sheet_visium.csv')
sample_sheet = pd.read_csv(sample_sheet_path)
sample_list = sample_sheet['section_id'].unique().tolist()
sample_id = 'P29_Vis'  # Example sample ID for testing