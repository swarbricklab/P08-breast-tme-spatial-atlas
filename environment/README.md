# Environment Documentation

This directory contains environment definition and snapshot files for reproducible analysis workflows.

## Environment Usage in Scripts

Each analysis script includes a line near the top specifying which environment was used. This helps ensure reproducibility and clarity when running or submitting scripts.

**Example:**
```python
# --- Define environment
# env: xen-hev
```

When submitting scripts as jobs (e.g., via HPC), the environment can also be set in the `params.sh` file:
```bash
# Environment to use
environment="xen-hev"
```
This allows job submission scripts to activate the correct environment automatically.

---

## Notes on Reproducibility

The analyses for Visium, Visium HD, and Xenium data were performed in separate environments (`spatial-pilot`, `vis-hd`, `xen-hev`) to ensure compatibility during development. **For reproducibility, you can use a single combined environment** (e.g., merging the dependencies from `spatial-pilot`, `vis-hd`, and `xen-hev`) to run all spatial transcriptomics analyses, though package versions may differ slightly from those used in the original analysis.

However, we recommend keeping the following environments separate as they are based on official package documentation and have specific version requirements:
- **cell2location**: For spatial deconvolution
- **tangram**: For mapping single-cell data to spatial coordinates
- **spatdat**: For SpatialData object creation and manipulation

Both Environment Definition Files (EDF) and Environment Snapshot Files (ESF) are provided for full reproducibility.

---

## Structure

- **EDF/**: Environment Definition Files. Use these YAML files to create new conda environments:
  ```bash
  conda env create -f environment/EDF/{env_file}.yaml
  ```
- **ESF/**: Environment Snapshot Files. These are exported from active environments to capture exact package versions:
  ```bash
  conda activate {env_name}
  conda env export --no-build > environment/ESF/{env_name}.yaml
  ```

## Environments

### Platform-Specific Analysis Environments (can be combined)

These environments were used separately during development but can be merged into a single environment for simplified reproduction:

#### spatial-pilot
- Used for Visium analysis
- Core packages: scanpy, squidpy, cell2location (inference only)

#### vis-hd
- Used for Visium HD analysis
- Core packages: scanpy, squidpy

#### xen-hev
- Used for Xenium analysis (e.g., TLS and HEV region identification)
- Core packages: scanpy, squidpy

### Specialized Tool Environments (recommended to keep separate)

These environments are based on official package documentation and should be kept separate due to specific version requirements:

#### cell2location
- Used for spatial deconvolution model training
- Based on official cell2location installation instructions

#### tangram
- Used for mapping single-cell reference data to spatial coordinates
- Based on official Tangram installation instructions

#### spatdat
- Used for creation and handling of SpatialData objects (Xenium data)
- Based on official SpatialData installation instructions

### Utility Environments

#### gseapy
- Used for gene set enrichment analysis with the gseapy package

#### shapely
- Used for scripts requiring spatial geometry operations, advanced annotation placement, and set visualization
- Key packages: `shapely`, `adjusttext`, `upsetplot`

---

For each environment, refer to the corresponding YAML file in EDF for dependencies and to ESF for a snapshot of the installed packages.
