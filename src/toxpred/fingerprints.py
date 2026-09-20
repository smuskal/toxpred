"""The one fingerprint this package uses for read-across.

Morgan, radius 2, 1,024 bits, which stands in for the Daylight fingerprint of
the original method and is open and standard.
"""
from __future__ import annotations

import numpy as np


def morgan(smiles, n_bits: int = 1024, radius: int = 2):
    """-> (bit matrix as bool, boolean mask of the SMILES that parsed)."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius,
                                                    fpSize=n_bits)
    rows, ok = [], []
    for s in smiles:
        m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
        if m is None:
            rows.append(np.zeros(n_bits, dtype=bool))
            ok.append(False)
            continue
        rows.append(np.asarray(gen.GetFingerprintAsNumPy(m), dtype=bool))
        ok.append(True)
    return np.stack(rows) if rows else np.zeros((0, n_bits), bool), \
        np.asarray(ok, dtype=bool)


def tanimoto(query_bits, ref_bits, block: int = 512):
    """Full query by reference Tanimoto matrix."""
    A = query_bits.astype(np.float32)
    B = ref_bits.astype(np.float32)
    na, nb = A.sum(1), B.sum(1)
    T = np.zeros((len(A), len(B)), dtype=np.float32)
    for s in range(0, len(A), block):
        e = min(s + block, len(A))
        inter = A[s:e] @ B.T
        union = na[s:e, None] + nb[None, :] - inter
        with np.errstate(invalid="ignore", divide="ignore"):
            T[s:e] = np.where(union > 0, inter / union, 0.0)
    return T
