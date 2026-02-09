import matplotlib.pyplot as plt
import numpy as np
from Bio import AlignIO

alignment = AlignIO.read("../data/aligned/spike_sequences.aln", "clustal")

RBD_START = 957
RBD_END = 1623

names = [rec.id for rec in alignment]

def percent_identity(s1, s2):
    m = v = 0
    for a,b in zip(s1, s2):
        if a == "-" or b == "-": continue
        v += 1
        if a == b: m += 1
    return 100*m/v

n = len(alignment)
matrix = np.zeros((n,n))

for i in range(n):
    for j in range(n):
        r1 = alignment[i].seq[RBD_START:RBD_END]
        r2 = alignment[j].seq[RBD_START:RBD_END]
        matrix[i,j] = percent_identity(r1, r2)

plt.figure(figsize=(6,5))
plt.imshow(matrix)
plt.xticks(range(n), names, rotation=45)
plt.yticks(range(n), names)
plt.colorbar(label="% identity")
plt.title("RBD divergence heatmap")
plt.tight_layout()
plt.savefig("../figures/rbd_heatmap1.png", dpi=300)
