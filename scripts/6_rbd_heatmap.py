import matplotlib.pyplot as plt
import numpy as np
from Bio import AlignIO
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent

alignment = AlignIO.read(BASE / "data/aligned/spike_proteins.aln", "clustal")

seqs = {}

for record in alignment:
    if record.id == "Wuhan":
        seqs["Wuhan"] = record.seq
    elif record.id == "Delta":
        seqs["Delta"] = record.seq
    elif record.id == "Omicron":
        seqs["Omicron"] = record.seq

# RBD region (protein positions)
start = 319 - 1
end = 541

positions = list(range(319, 541))
variants = ["Delta", "Omicron"]

matrix = []

for var in variants:
    row = []
    for i in range(start, end):
        w = seqs["Wuhan"][i]
        v = seqs[var][i]

        if w != v and w not in "-X" and v not in "-X":
            row.append(1)
        else:
            row.append(0)

    matrix.append(row)

matrix = np.array(matrix)

plt.figure(figsize=(16,4))

im = plt.imshow(
    matrix,
    aspect="auto",
    interpolation="nearest",
    cmap="Reds"   # cleaner for mutation presence
)

# Axis labels
plt.yticks(range(len(variants)), variants, fontsize=12)
plt.xticks(
    range(0, len(positions), 15),
    positions[::15],
    fontsize=11,
    rotation=0
)

plt.xlabel("RBD amino acid position", fontsize=13)
plt.ylabel("Variant", fontsize=13)
plt.title("RBD mutation landscape relative to Wuhan", fontsize=15, pad=10)

# Subtle grid for readability
plt.grid(axis="x", linestyle="--", alpha=0.25)

# Horizontal separator between variants
plt.axhline(y=0.5, color="black", linewidth=1.5)

# Colorbar (clean)
cbar = plt.colorbar(im, fraction=0.02, pad=0.02)
cbar.set_label("Mutation present", fontsize=12)

plt.tight_layout()

fig_path = BASE / "figures/rbd_heatmap.png"
fig_path.parent.mkdir(exist_ok=True)

plt.savefig(fig_path, dpi=400)
plt.show()

print("Saved → rbd_heatmap.png")
