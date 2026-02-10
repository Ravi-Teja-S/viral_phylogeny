from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from pathlib import Path
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent.parent

alignment = AlignIO.read(
    BASE / "data/aligned/spike_proteins.aln",
    "fasta"
)

calculator = DistanceCalculator("blosum62")
distance_matrix = calculator.get_distance(alignment)

constructor = DistanceTreeConstructor()
tree = constructor.nj(distance_matrix)

tree.root_with_outgroup("Bat_outgroup")
tree.ladderize()

tree_file = BASE / "results/phylogeny_tree.nwk"
tree_file.parent.mkdir(exist_ok=True)
Phylo.write(tree, tree_file, "newick")

fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot(1,1,1)
Phylo.draw(tree, axes=ax, do_show=False)

plt.tight_layout()
plt.savefig(BASE / "figures/phylogeny.png", dpi=300)
plt.show()
