"""The consortium method: a compound inherits from the compounds it resembles.

A query is compared to every compound in a reference set. Those at or above a
similarity cutoff form its voting consortium. A prediction is issued only when
the consortium has at least a minimum number of members, and is the unweighted
mean of their values for a continuous property or the majority of their labels
for a categorical one. Weighting votes by similarity buys nothing and is not
done.

WHEN THERE IS NO CLOSE ENOUGH NEIGHBOR THE METHOD DECLINES. That is the point.
Accuracy and coverage are separate quantities and both are reported, always.
Raising the cutoff raises accuracy and lowers the share of compounds you get an
answer for, and there is no setting that escapes the trade.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .fingerprints import morgan, tanimoto

# GHS categories, which are assigned in the source data rather than by us.
#   RED    at or under 50 mg/kg
#   AMBER  50 to 2000
#   GREEN  over 2000
GHS_LIGHT = {1: "RED", 2: "RED", 3: "AMBER", 4: "AMBER", 5: "GREEN"}


@dataclass
class Prediction:
    query: str
    call: object = None
    confidence: float = float("nan")
    members: int = 0
    nearest: float = float("nan")
    answered: bool = False
    exact_match: object = None
    neighbors: list = field(default_factory=list)

    def __str__(self):
        if not self.answered:
            return ("%s: no call, nearest reference compound at %.3f"
                    % (self.query, self.nearest))
        extra = ""
        if self.exact_match is not None:
            extra = "; an identical structure in the reference set says %s" \
                % self.exact_match
        return ("%s: %s, %d voters agreeing %.0f%%, nearest %.3f%s"
                % (self.query, self.call, self.members,
                   100 * self.confidence, self.nearest, extra))


class Consortium:
    """A reference set plus the voting rule.

    reference : DataFrame with a `smiles` column and one label or value column
    value_col : the column to inherit
    kind      : "label" for a majority vote, "value" for an unweighted mean
    """

    def __init__(self, reference: pd.DataFrame, value_col: str,
                 kind: str = "label", n_bits: int = 1024):
        if "smiles" not in reference.columns:
            raise ValueError("reference needs a 'smiles' column")
        if value_col not in reference.columns:
            raise ValueError("reference has no column %r" % value_col)
        bits, ok = morgan(reference.smiles.tolist(), n_bits=n_bits)
        self.reference = reference[ok].reset_index(drop=True)
        self.bits = bits[ok]
        self.values = self.reference[value_col].values
        self.kind = kind
        self.n_bits = n_bits

    def __len__(self):
        return len(self.reference)

    @property
    def fingerprint(self) -> str:
        """A short hash of the reference set actually in use.

        Quote it beside any number this produces. Two runs that agree on this
        string were answering from the same reference set; two that do not were
        not, however similar the rest of the setup looked.
        """
        import hashlib
        h = hashlib.sha256()
        for smi, val in zip(self.reference.smiles, self.values):
            h.update(("%s\t%s\n" % (smi, val)).encode())
        return h.hexdigest()[:12]

    def predict(self, smiles, cutoff: float = 0.5, min_members: int = 1,
                keep_neighbors: int = 0):
        """-> list of Prediction, one per query, in order."""
        smiles = [smiles] if isinstance(smiles, str) else list(smiles)
        qbits, ok = morgan(smiles, n_bits=self.n_bits)
        T = tanimoto(qbits, self.bits)
        out = []
        for i, s in enumerate(smiles):
            p = Prediction(query=s)
            if not ok[i]:
                out.append(p)
                continue
            sims = T[i]
            p.nearest = float(sims.max()) if len(sims) else float("nan")
            m = sims >= cutoff
            p.members = int(m.sum())
            if p.members < min_members:
                out.append(p)
                continue
            vals = self.values[m]
            if self.kind == "value":
                p.call = float(np.mean(vals))
                spread = float(np.std(vals))
                p.confidence = float(1.0 / (1.0 + spread))
            else:
                names, counts = np.unique(vals, return_counts=True)
                order = np.argsort(-counts)
                # a tie is not a majority; the method declines rather than
                # breaking it on whichever label numpy happened to sort first
                if len(order) > 1 and counts[order[0]] == counts[order[1]]:
                    p.confidence = float(counts[order[0]] / counts.sum())
                    out.append(p)
                    continue
                j = int(order[0])
                p.call = names[j]
                p.confidence = float(counts[j] / counts.sum())
            p.answered = True
            if p.nearest >= 0.999:
                k = int(np.argmax(sims))
                p.exact_match = self.values[k]
            if keep_neighbors:
                order = np.argsort(-sims)[:keep_neighbors]
                p.neighbors = [(self.reference.smiles.iloc[int(k)],
                                float(sims[int(k)]), self.values[int(k)])
                               for k in order if sims[int(k)] >= cutoff]
            out.append(p)
        return out


class TrafficLight(Consortium):
    """Red, amber and green for rat oral acute toxicity, from CATMoS."""

    @classmethod
    def from_catmos(cls, catmos_csv, **kw):
        d = pd.read_csv(catmos_csv)
        d["light"] = pd.to_numeric(d.ghs_cat, errors="coerce").map(GHS_LIGHT)
        d = d[d.light.notna() & d.smiles.notna()].reset_index(drop=True)
        return cls(d, "light", kind="label", **kw)


class Endpoint(Consortium):
    """One of the standardized binary endpoints."""

    @classmethod
    def from_toxric(cls, csv_path, **kw):
        d = pd.read_csv(csv_path)
        d = d.rename(columns={"Canonical SMILES": "smiles",
                              "Toxicity Value": "value"})
        d = d[d.smiles.notna()].drop_duplicates("smiles").reset_index(drop=True)
        d["value"] = pd.to_numeric(d.value, errors="coerce")
        d = d[np.isfinite(d.value)].reset_index(drop=True)
        d["call"] = np.where(d.value == 1, "TOXIC", "clear")
        return cls(d, "call", kind="label", **kw)
