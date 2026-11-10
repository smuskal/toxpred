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


def largest_fragment(smiles):
    """The largest fragment of a structure, canonical, or None if unreadable.

    ONE RULE, and it is the same one `pharmcast.features._largest` applies:
    max on (heavy atom count, canonical SMILES). The SMILES is the tie-break
    and it is not decoration -- `max` returns the FIRST maximum, so heavy-atom
    count alone makes the answer depend on the order the fragments were typed.
    Measured: "CCO.CCN" keeps CCO and "CCN.CCO" keeps CCN. Breaking the tie on
    the canonical SMILES makes the choice depend only on the chemistry, so two
    spellings of one 1:1 co-former salt score the same molecule.

    This is a second implementation of a rule that also lives in
    pharmcast.features, which is a real risk, so it is GATED rather than
    trusted: selftest.salt_handling() checks the two agree on every
    multi-fragment structure in the reference sets, and the deploy gate runs
    it. They were verified identical on all 423 multi-fragment TOXRIC rows and
    on CATMoS (which has none) on 9 Nov 2026. The copy exists because the
    toxpred package is installable without pharmcast.
    """
    from rdkit import Chem
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    if "." not in smiles:
        return Chem.MolToSmiles(mol)
    try:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    except Exception:
        return None
    if len(frags) <= 1:
        return Chem.MolToSmiles(mol)
    return Chem.MolToSmiles(max(frags, key=lambda f: (f.GetNumHeavyAtoms(),
                                                      Chem.MolToSmiles(f))))


def desalt_and_resolve(smiles, values):
    """-> (smiles, values, report) with both sides reduced to one structure.

    THE ORDER IS THE POINT: strip FIRST, then group. De-duplicating on the raw
    SMILES, which is what this used to do, left a salt and its own parent
    standing as two rows, so one compound voted twice and inflated the
    consortium it was a member of.

    A group whose labels all agree collapses to one row. A group whose labels
    DISAGREE is dropped whole: after stripping, those rows are one structure
    carrying two opposite measurements, there is no majority to take, and a
    reference set should not assert what its own evidence contradicts.
    `Consortium.predict` declines on a tie for the same reason. Letting the
    parent row win, or taking the majority, were both measured on 9 Nov 2026:
    the median difference across the 28 endpoints is 0.0000 either way, but
    under `parent` SR-ARE's operating cutoff moves, which is a selection
    artifact rather than a result. Dropping is the only policy that leaves
    every operating cutoff where it was.

    The report is returned rather than logged so a caller can state the cost
    instead of the cost being invisible.
    """
    groups = {}
    unreadable = 0
    for smi, val in zip(smiles, values):
        key = largest_fragment(smi)
        if key is None:
            unreadable += 1
            continue
        groups.setdefault(key, []).append(val)
    out_smi, out_val = [], []
    merged = conflicts = conflict_rows = 0
    for key, vals in groups.items():
        distinct = set(vals)
        if len(distinct) == 1:
            if len(vals) > 1:
                merged += len(vals) - 1
            out_smi.append(key)
            out_val.append(vals[0])
        else:
            conflicts += 1
            conflict_rows += len(vals)
    report = dict(rows_in=len(smiles), rows_out=len(out_smi),
                  rows_merged_away=merged, conflict_groups_dropped=conflicts,
                  rows_lost_to_conflicts=conflict_rows,
                  unreadable=unreadable)
    return out_smi, out_val, report

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
        # The secondary sort key for the neighbor list, built once here rather
        # than per query. See the lexsort in predict(): it is what makes two
        # neighbors at the same similarity come out in the same order every
        # time, in this process and in the next one.
        self._nb_key = np.asarray(self.reference.smiles.tolist(), dtype=object)

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
                # A TIE IS BROKEN ON THE SMILES, SO THE LIST IS REPRODUCIBLE.
                # `np.argsort(-sims)` is quicksort, which is not stable, so two
                # neighbors at the same similarity came out in whichever order
                # the partition happened to leave them. Measured 7 Nov 2026:
                # imatinib's list swapped tied pairs between two revisions of
                # this service while every call, confidence and voter count was
                # identical. Nothing was wrong with the answer, but a
                # before-and-after diff could never be read, because noise and
                # a real change looked the same.
                #
                # lexsort takes the LAST key as primary, so this is similarity
                # descending, then canonical SMILES ascending. That is the same
                # tie-break as features._largest(), which is deliberate: one
                # convention for every ordering decision in this stack.
                order = np.lexsort((self._nb_key, -sims))[:keep_neighbors]
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
        # DESALTED, AND DELIBERATELY NOT DE-DUPLICATED.
        #
        # Every SMILES is reduced to its largest fragment so the reference set
        # matches the query side. That is a no-op on the data that ships:
        # CATMoS carries ZERO multi-fragment rows, measured 9 Nov 2026, so the
        # traffic light is bit-for-bit what it was. It is here so a future
        # CATMoS release carrying a salt cannot quietly reintroduce the
        # asymmetry this stack just removed.
        #
        # It does NOT go through desalt_and_resolve, and that is the point.
        # That helper also merges duplicate structures and drops contradictory
        # groups, which this set has never had done to it: running it here
        # took the acute traffic light from 11,472 rows to 10,628, merging 591
        # and dropping 97 contradictory groups, and moved the reference
        # fingerprint. That is a change to the headline acute deliverable,
        # it is not desalting, and it was not what anyone asked for. Caught
        # before it shipped on 9 Nov 2026. De-duplicating CATMoS may well be
        # worth doing; it is its own decision, with its own measurement.
        # A SINGLE-FRAGMENT ROW IS PASSED THROUGH VERBATIM, not re-canonicalized.
        # largest_fragment() round-trips every structure through MolToSmiles,
        # which rewrites the spelling even when there was only ever one
        # fragment. Morgan bits are invariant to spelling, so no call or
        # confidence moves either way -- but the reference `fingerprint` hashes
        # the SMILES strings, and the neighbor lists shown to a reader print
        # them. Re-canonicalizing all 11,472 rows moved the fingerprint from
        # 8216f99c4ee5 to 3b920e555058 while every prediction stayed the same,
        # which is the worst kind of diff: a changed identity with no changed
        # behavior behind it. Touching only the rows that really carry a second
        # fragment makes this a true no-op on the data that ships.
        kept, dropped = [], 0
        for smi, light in zip(d.smiles.tolist(), d.light.tolist()):
            if isinstance(smi, str) and "." not in smi:
                kept.append((smi, light))
                continue
            one = largest_fragment(smi)
            if one is None:
                dropped += 1
                continue
            kept.append((one, light))
        d = pd.DataFrame(kept, columns=["smiles", "light"])
        obj = cls(d, "light", kind="label", **kw)
        obj.desalt_report = dict(rows_in=len(kept) + dropped,
                                 rows_out=len(kept), unreadable=dropped,
                                 deduplicated=False)
        return obj


class Endpoint(Consortium):
    """One of the standardized binary endpoints."""

    @classmethod
    def from_toxric(cls, csv_path, **kw):
        d = pd.read_csv(csv_path)
        d = d.rename(columns={"Canonical SMILES": "smiles",
                              "Toxicity Value": "value"})
        # drop_duplicates on the RAW smiles still runs, and it still only
        # removes rows that were already identical as shipped. The real
        # de-duplication now happens on the STRIPPED structure below, which is
        # what catches a salt sitting beside its own parent.
        d = d[d.smiles.notna()].drop_duplicates("smiles").reset_index(drop=True)
        d["value"] = pd.to_numeric(d.value, errors="coerce")
        d = d[np.isfinite(d.value)].reset_index(drop=True)
        d["call"] = np.where(d.value == 1, "TOXIC", "clear")
        # BOTH SIDES ARE DESALTED FROM 9 Nov 2026. Until then these sets were
        # used as TOXRIC ships them -- 2,236 of 168,075 rows carrying a
        # counterion -- and the query was not reduced either, because reducing
        # only the query would not have matched them. That cost answers:
        # propranolol called AMBER on 15 endpoints while its hydrochloride
        # declined on 14, the counterion having dropped even its own parent in
        # the reference set from 1.000 to 0.971 and pushed the rest below the
        # voting cutoff.
        smis, vals, rep = desalt_and_resolve(d.smiles.tolist(),
                                             d.call.tolist())
        d = pd.DataFrame({"smiles": smis, "call": vals})
        obj = cls(d, "call", kind="label", **kw)
        obj.desalt_report = rep
        return obj
