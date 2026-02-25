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

### cell2location
- Used for deconvolution of Visium data with cell2location package.

### gseapy
- Used for gene set enrichment analysis with the gseapy package.

### spatdat
- Used for creation and handling of SpatialData objects.

### spatial-pilot
- Used for Visium analysis.

### vis-hd
- Used for Visium HD analysis.

### xen-hev
- Used for Xenium analysis (e.g., TLS and HEV region identification).

### shapely
- Used for scripts that require spatial operations with the `shapely` package, advanced annotation placement with `adjusttext`, and set visualization with `upsetplot`.
- Key packages: `shapely`, `adjusttext`, `upsetplot`

---

For each environment, refer to the corresponding YAML file in EDF for dependencies and to ESF for a snapshot of the installed packages.
