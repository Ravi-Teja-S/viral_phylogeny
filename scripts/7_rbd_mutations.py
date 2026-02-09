from Bio import AlignIO

alignment = AlignIO.read("../data/aligned/spike_proteins.aln", "clustal")

# identify sequences
wuhan = None
omicron = None

for record in alignment:
    if "NC_045512" in record.id or "Wuhan" in record.id:
        wuhan = record.seq
    if "Omicron" in record.id or "OM" in record.id:
        omicron = record.seq

# RBD region in spike protein (approx 319–541)
rbd_start = 319 - 1
rbd_end = 541

print("\nOmicron RBD mutations vs Wuhan\n")
print("-"*50)
print("Pos  Wuhan  Omicron  Mutation")
print("-"*50)

pos = 319
for i in range(rbd_start, rbd_end):
    w = wuhan[i]
    o = omicron[i]

    if w != o and w not in "-X" and o not in "-X":
        print(f"{pos:<4} {w:<6} {o:<8} {w}{pos}{o}")

    pos += 1
