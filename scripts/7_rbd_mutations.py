from Bio import AlignIO

alignment = AlignIO.read(
    "../data/aligned/spike_proteins.aln",
    "fasta"
)

# ---- Load all sequences into dict ----
seqs = {rec.id: rec.seq for rec in alignment}

# ---- Find Wuhan reference robustly ----
ref_name = None
for name in seqs:
    if "Wuhan" in name or "NC_045512" in name:
        ref_name = name
        break

if ref_name is None:
    raise ValueError("Wuhan reference not found in alignment")

ref = seqs[ref_name]

# ---- RBD region (protein coords) ----
rbd_start = 319 - 1
rbd_end = 541

print(f"\nReference sequence → {ref_name}")
print("=" * 60)

# ---- Compare all variants to Wuhan ----
for name, seq in seqs.items():

    if name == ref_name:
        continue

    print(f"\n{name} RBD mutations vs Wuhan")
    print("-" * 55)
    print("Pos  Wuhan  Variant  Mutation")
    print("-" * 55)

    pos = 319
    mutation_count = 0

    for i in range(rbd_start, rbd_end):
        w = ref[i]
        v = seq[i]

        if w != v and w not in "-X" and v not in "-X":
            print(f"{pos:<4} {w:<6} {v:<8} {w}{pos}{v}")
            mutation_count += 1

        pos += 1

    if mutation_count == 0:
        print("No RBD mutations detected")

print("\nDone.")
