"""The grouping used to compare donors: one column of the cell table, an optional reference group, and the order of the groups.

A column can group donors when every donor has one value in it, with 2 to MAX_LEVELS values overall. The choice is saved in the project's cache
(groups.json); until the user chooses, a column is suggested from its name. A reference is suggested only when a group is plainly named as one
(control, healthy ...); otherwise there is none and no "versus reference" test is run. A trend test across the groups needs `ordered`, which only
the user can set, because the app cannot know that the groups form a scale.
"""
import json, re
import numpy as np, pandas as pd

MAX_LEVELS = 30
NO_GROUP = "All donors"                      # the one group when no column is chosen
NAME_HINTS = ("subtype", "group", "condition", "diagnosis", "disease", "status", "treatment", "genotype")      # likeliest first
REFERENCE_NAME = re.compile(r"controls?|ctrl|healthy|normal|reference|untreated|vehicle|wild[ _-]?type|wt", re.I)
ID_NAME = re.compile(r"(^|[._-])(id|barcode|uuid)([._-]|$)|blinded|merged_from", re.I)


def natural(s):
    """Sort key: "Group 10" after "Group 2"."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(s))]


def candidates(store):
    """{column: {sample: value}} of the columns that can group this project's donors (none with a single donor)."""
    df, samples = store.df, store.samples
    donors = pd.factorize(pd.Series([store.donor[s] for s in samples]))[0]
    n_donors = int(donors.max()) + 1
    if n_donors < 2:
        return {}
    starts = np.array([store.slices[s][0] for s in samples])
    used = {c for c in store.map.values() if c}
    out = {}
    for col in df.columns:
        if col in used or df[col].dtype.kind == "f":
            continue
        codes, values = pd.factorize(df[col])
        if not 2 <= len(values) <= MAX_LEVELS:
            continue
        lo, hi = np.minimum.reduceat(codes, starts), np.maximum.reduceat(codes, starts)
        if (lo < 0).any() or (lo != hi).any():       # missing, or more than one value inside a tissue
            continue
        if (pd.Series(lo).groupby(donors).nunique() > 1).any():
            continue
        levels = np.unique(lo)
        if len(levels) < 2 or (len(levels) == n_donors and ID_NAME.search(str(col))):      # an id of the donor, not a grouping
            continue
        out[str(col)] = {s: str(values[c]) for s, c in zip(samples, lo)}
    return out


class Groups:
    def __init__(self, store):
        self.options = candidates(store)
        self.file = store.cache_dir / "groups.json"
        self.column, self.reference, self.order, self.ordered = None, None, [], False
        saved = None
        if self.file.is_file():
            try:
                saved = json.loads(self.file.read_text())
            except ValueError:
                pass
        if isinstance(saved, dict) and (saved.get("column") is None or saved["column"] in self.options):
            self._apply(saved.get("column"), saved.get("reference"), saved.get("order"), bool(saved.get("ordered")))
        else:
            self._apply(self._suggest(), "auto", None, False)

    def _suggest(self):
        for hint in NAME_HINTS:
            for c in self.options:
                if hint in c.lower():
                    return c
        return None

    def levels(self, column):
        return sorted(set(self.options[column].values()), key=natural)

    def _apply(self, column, reference, order, ordered):
        """Set the grouping; a saved reference or order that no longer fits the column falls back to the default."""
        if column is None:
            self.column, self.reference, self.order, self.ordered = None, None, [], False
            return
        lv = self.levels(column)
        if reference == "auto":
            reference = next((x for x in lv if REFERENCE_NAME.fullmatch(x.strip())), None)
        self.column, self.reference = column, reference if reference in lv else None
        if order is None or sorted(order) != sorted(lv):
            order = ([self.reference] if self.reference else []) + [x for x in lv if x != self.reference]
        self.order, self.ordered = list(order), bool(ordered) and len(lv) >= 3

    def set(self, column, reference=None, order=None, ordered=False):
        """The user's choice; ValueError if it does not fit this project. Saved in the cache when the cache can be written."""
        if column is not None:
            if column not in self.options:
                raise ValueError(f"{column!r} cannot group the donors")
            lv = self.levels(column)
            if reference is not None and reference != "auto" and reference not in lv:
                raise ValueError(f"{reference!r} is not a group of {column}")
            if order is not None and sorted(order) != sorted(lv):
                raise ValueError(f"the order must list each group of {column} once")
        self._apply(column, reference, order, ordered)
        try:
            self.file.write_text(json.dumps({"column": self.column, "reference": self.reference, "order": self.order, "ordered": self.ordered}))
        except OSError:
            pass

    def of(self, sample):
        return self.options[self.column][sample] if self.column else NO_GROUP

    def describe(self):
        return {"column": self.column, "reference": self.reference, "order": self.order, "ordered": self.ordered,
                "options": [{"column": c, "levels": self.levels(c)} for c in self.options]}
