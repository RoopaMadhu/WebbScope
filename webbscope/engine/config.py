"""Column mapping: which metadata columns play which role. Overridable via config.json."""
import json
from pathlib import Path

# Candidate column names used to suggest a mapping when none is configured.
CANDIDATES = {
    "cell_id": ["cell_id", "cell", "barcode", "CellID"],
    "x": ["x", "x_centroid", "X", "center_x"],
    "y": ["y", "y_centroid", "Y", "center_y"],
    "sample": ["sample_id_blinded", "sample_id", "sample", "tissue_id"],
    "donor": ["donor_id_blinded", "donor_id", "donor"],
    "cell_type": ["celltype.coarse.lvl2", "celltype.coarse", "cell_type", "celltype", "cell_type_coarse"],
    "cell_state": ["cell_type_fine_lvl1", "cell_state", "cell_type_fine", "state"],
    "region_lvl1": ["annotation_broad", "region_lvl1", "region1"],
    "region_lvl2": ["annotation_clean", "region_lvl2", "region2"],
    "region_lvl3": ["annotation_fine_lvl1", "region_lvl3", "region3", "tessera_niche"],
    "lib_size": ["lib_size"],
}

REGION_ROLES = ["region_lvl1", "region_lvl2", "region_lvl3"]
CATEGORICAL_ROLES = ["cell_type", "cell_state"] + REGION_ROLES


def find_columns(columns, word):
    """Columns whose name contains `word` (any case), the likeliest id column first: names that look like a type / count / number of something
    last, then names with "id" or "name" before others, then the shortest name, then table order."""
    hit = [(i, str(c)) for i, c in enumerate(columns) if word in str(c).lower()]
    hit.sort(key=lambda t: (any(w in t[1].lower() for w in ("type", "count", "n_", "num")), not any(w in t[1].lower() for w in ("id", "name")), len(t[1]), t[0]))
    return [c for _, c in hit]


def suggest(columns):
    cols = set(columns)
    out = {role: next((c for c in cands if c in cols), None) for role, cands in CANDIDATES.items()}
    for role in ("sample", "donor"):                      # no usual name present: any column that has the word in its name
        if out[role] is None:
            f = find_columns(columns, role)
            out[role] = f[0] if f else None
    return out


def load_mapping(columns, path=None):
    mapping = suggest(columns)
    if path and Path(path).exists():
        mapping.update({k: v for k, v in json.loads(Path(path).read_text()).items() if v is None or v in columns})
    return mapping
