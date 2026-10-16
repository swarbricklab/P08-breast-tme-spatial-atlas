# --- Discrete-TLS definition (Extended Data Fig. 3b)
# Derives "discrete TLS" from the permissive per-spot TLS grouping produced by
# 02_group_TLS_high (contiguous groups of >=2 directly adjacent TLS-high spots,
# score > 0.01, in the 'group_singleton_removed' column).
#
# Stringent Extended Data Fig. 3b definition (per-sample quantification + spatial analysis):
#   1. Keep only groups of >= min_spots directly adjacent TLS-high spots.
#   2. Merge surviving groups separated by a single spot (nearest centre-to-centre
#      spot distance <= ~2 spot pitches).

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist


def _estimate_pitch_um(coords_arr):
    """Median nearest-neighbour spot spacing (~100 µm for Visium)."""
    if len(coords_arr) < 2:
        return 100.0
    tree = cKDTree(coords_arr)
    dd, _ = tree.query(coords_arr, k=2)
    pitch = float(np.median(dd[:, 1]))
    return pitch if np.isfinite(pitch) and pitch > 0 else 100.0


def define_discrete_tls(groups_df, coords_um, group_col='group_singleton_removed',
                        min_spots=3, merge_gap_pitches=2.2):
    """Return a per-spot Series (index = barcode) of discrete-TLS labels.

    - groups_df: per-spot table with `group_col` (permissive grouping; NaN where
      the spot is not in a >=2-spot group).
    - coords_um: per-spot centre coordinates in µm (index = barcode), covering at
      least all TLS-high spots; used to estimate the spot pitch and to measure
      inter-group distances.
    - min_spots: minimum directly-adjacent spots for a group to survive.
    - merge_gap_pitches: merge two surviving groups if their nearest centre-to-centre
      spot distance is <= merge_gap_pitches * pitch (a single intervening spot is
      ~2 pitches; 2.2 admits both the colinear (2 pitch) and diagonal (sqrt(3) pitch)
      one-spot gaps while excluding a two-spot gap).

    Spots not in a surviving/merged discrete TLS are omitted from the result.
    Discrete-TLS ids are 'TLS1', 'TLS2', ... ordered by size (largest first).
    """
    labels = groups_df[group_col].dropna()
    if labels.empty:
        return pd.Series(dtype=object)

    sizes = labels.value_counts()
    survivor_ids = sizes.index[sizes >= min_spots]
    surv = labels[labels.isin(survivor_ids)]
    group_ids = list(surv.unique())
    n = len(group_ids)
    if n == 0:
        return pd.Series(dtype=object)

    pitch = _estimate_pitch_um(coords_um.values)
    thresh = merge_gap_pitches * pitch

    coords_by_group = {
        g: coords_um.loc[surv.index[surv == g]].values for g in group_ids
    }

    # Union-find merge of surviving groups within `thresh`
    parent = {g: g for g in group_ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(n):
        gi = group_ids[i]
        ci = coords_by_group[gi]
        for j in range(i + 1, n):
            gj = group_ids[j]
            if find(gi) == find(gj):
                continue
            if cdist(ci, coords_by_group[gj]).min() <= thresh:
                union(gi, gj)

    roots = {g: find(g) for g in group_ids}
    comp_size = {}
    for g in group_ids:
        comp_size[roots[g]] = comp_size.get(roots[g], 0) + len(coords_by_group[g])
    root_to_id = {
        root: f'TLS{k + 1}'
        for k, root in enumerate(sorted(comp_size, key=lambda r: -comp_size[r]))
    }

    return surv.map(lambda g: root_to_id[roots[g]])
