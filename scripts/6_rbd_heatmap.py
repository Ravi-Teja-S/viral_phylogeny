import matplotlib.pyplot as plt
import numpy as np
from Bio import AlignIO
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent

alignment = AlignIO.read(
    BASE / "data/aligned/spike_proteins.aln",
    "fasta"
)

# ---- Pick reference dynamically (Wuhan-like) ----
REFERENCE_NAME = "Wuhan"

seqs = {rec.id: rec.seq for rec in alignment}

if REFERENCE_NAME not in seqs:
    raise ValueError(f"Reference '{REFERENCE_NAME}' not found in alignment")

ref = seqs[REFERENCE_NAME]

# ---- All other variants automatically ----
variants = [k for k in seqs if k != REFERENCE_NAME]

# ---- RBD region (protein coordinates) ----
start = 319 - 1
end = 541

positions = list(range(319, 541))

matrix = []

for var in variants:
    row = []
    for i in range(start, end):
        r = ref[i]
        v = seqs[var][i]

        if r != v and r not in "-X" and v not in "-X":
            row.append(1)
        else:
            row.append(0)
    matrix.append(row)

matrix = np.array(matrix)

# ---- Plot ----
plt.figure(figsize=(18, max(2, len(variants) * 0.6)))

im = plt.imshow(
    matrix,
    aspect="auto",
    interpolation="nearest",
    cmap="Reds"
)

plt.yticks(range(len(variants)), variants, fontsize=11)

step = max(1, len(positions) // 15)
plt.xticks(
    range(0, len(positions), step),
    positions[::step],
    fontsize=10
)

plt.xlabel("RBD amino acid position", fontsize=13)
plt.ylabel("Variants", fontsize=13)
plt.title("RBD mutation landscape relative to Wuhan", fontsize=15)

plt.grid(axis="x", linestyle="--", alpha=0.25)

cbar = plt.colorbar(im, fraction=0.025, pad=0.02)
cbar.set_label("Mutation present", fontsize=12)

plt.tight_layout()

fig_path = BASE / "figures/rbd_heatmap.png"
fig_path.parent.mkdir(exist_ok=True)

plt.savefig(fig_path, dpi=400)
plt.show()

print("Saved →", fig_path)
