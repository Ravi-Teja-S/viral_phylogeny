import math

# Transition pairs (RNA/DNA)
TRANSITIONS = {
    ("A", "G"), ("G", "A"),
    ("C", "T"), ("T", "C")
}

def kimura_distance(seq1, seq2):
    transitions = 0
    transversions = 0
    valid_sites = 0

    for a, b in zip(seq1.upper(), seq2.upper()):
        if a == "-" or b == "-":
            continue  # skip gaps

        if a == b:
            valid_sites += 1
            continue

        valid_sites += 1

        if (a, b) in TRANSITIONS:
            transitions += 1
        else:
            transversions += 1

    if valid_sites == 0:
        return 0

    P = transitions / valid_sites
    Q = transversions / valid_sites

    try:
        d = -0.5 * math.log(1 - 2*P - Q) - 0.25 * math.log(1 - 2*Q)
        return d
    except ValueError:
        return float("inf")   # handles extreme divergence
from Bio.Phylo.TreeConstruction import DistanceMatrix

def kimura_matrix(alignment):
    names = [rec.id for rec in alignment]
    matrix = []

    for i in range(len(alignment)):
        row = []
        for j in range(i + 1):
            d = kimura_distance(alignment[i].seq, alignment[j].seq)
            row.append(d)
        matrix.append(row)

    return DistanceMatrix(names, matrix)
