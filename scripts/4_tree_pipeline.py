import matplotlib.pyplot as plt
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from kimura_dist import kimura_matrix
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent


# ---------- Load alignment ----------
alignment = AlignIO.read(str(BASE / "data/aligned/spike_sequences.aln"), "clustal")
print("Alignment loaded with", len(alignment), "sequences")


# ---------- Distance calculation ----------
# Kimura 2-parameter model (good for nucleotide evolution)
distance_matrix = kimura_matrix(alignment)

print("\nDistance Matrix:\n")
print(distance_matrix)


# ---------- Neighbor Joining tree ----------
constructor = DistanceTreeConstructor()
tree = constructor.nj(distance_matrix)


# ---------- Root using bat coronavirus ----------
tree.root_with_outgroup("Bat_outgroup")
tree.ladderize()


#----------- Save tree file( Newick format) ---------
tree_file = BASE / "results/phylogeny_tree.nwk"
tree_file.parent.mkdir(exist_ok=True)

Phylo.write(tree,tree_file,"newick")
print("Tree file saved -> ",tree_file)

# ---------- Visualization ----------
fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot(1,1,1)

Phylo.draw(tree, axes=ax, do_show=False)

plt.title("Spike gene phylogeny of SARS-CoV-2 variants (NJ + Kimura)")
plt.xlabel("Genetic distance (substitutions per site)")
plt.tight_layout()

fig_path = BASE / "figures/phylogeny.png"
fig_path.parent.mkdir(exist_ok=True)

plt.savefig(fig_path, dpi=300)

print("\nSaved tree → phylogeny.png")
