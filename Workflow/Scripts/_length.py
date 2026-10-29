"""Gene-length correction shared by both arms' tidy steps.

Longer genes score better by chance (each pair's top 1% picks genes 1.39x its median
length). Each pair's energies are ranked within a length bin and mapped onto the middle
bin's distribution, so the value stays in kcal/mol; measured flatness 0.989 against 0.952
for the best fitted form, which the energy gate's truncation biases. An added ranking
column, never a gate.
"""

import numpy as np
import pandas as pd

# Length bins per pair. Flatness is 0.978-1.003 anywhere from 5 to 50 bins; ten keeps ~230
# genes per bin at this project's panel sizes.
LENGTH_BINS = 10

# Below this many hits a pair keeps its uncorrected value, so the column never holds NaN.
MIN_HITS = LENGTH_BINS * 5

# One tidy output is one (sample x variant), so grouping by miRNA is grouping by pair.
GROUP_COL = "miRNA"
LENGTH_COL = "Gene_length"


#----- Quantile-normalizes one pair's energies to its middle length bin -----#
def _correct_group(energy, length):
    if len(energy) < MIN_HITS:
        return energy

    binned = pd.qcut(length, LENGTH_BINS, labels=False, duplicates="drop")
    reference = energy[binned == binned.max() // 2]
    if len(reference) < 5:
        return energy

    # Rank inside the hit's own length bin, then read that rank off the reference bin.
    # Monotone within any single gene (one length, one bin), so the per-pair cap and
    # intersect.py's per-(miRNA, Gene) best-hit choice are provably unaffected.
    # np.quantile takes the reference unsorted — it computes order statistics itself.
    percentile = pd.Series(energy).groupby(binned).rank(pct=True, method="average")
    return np.quantile(reference, percentile.to_numpy())


#----- The added column's name. Spelled once because both tidy steps' OUTPUT_COLUMNS and
#      intersect's CONSENSUS_COLUMNS select on it, and each of those drops a name it does not
#      recognise silently (filtered out, or reindexed to all-NA) -----#
def corrected_name(energy_col):
    return f"{energy_col}_corrected"


#----- Adds `<energy_col>_corrected`, computed independently within each miRNA -----#
def add_length_corrected(df, energy_col):
    energy = df[energy_col].to_numpy(dtype=float)
    length = df[LENGTH_COL].to_numpy()
    corrected = np.full(len(df), np.nan)

    # Positional scatter via .indices: correct on tidy_intarna's sort-scrambled index, and
    # unlike groupby.apply it never reshapes a same-length result into a frame.
    for _, index in df.groupby(GROUP_COL, sort=False).indices.items():
        corrected[index] = _correct_group(energy[index], length[index])

    return df.assign(**{corrected_name(energy_col): corrected})
