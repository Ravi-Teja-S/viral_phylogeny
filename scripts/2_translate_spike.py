from Bio import SeqIO

with open("../data/proteins/spike_proteins.fasta", "w") as out:
    for record in SeqIO.parse("../data/raw/spike_sequences.fasta", "fasta"):
        protein = record.seq.translate()
        protein = protein.rstrip("*")

        out.write(f">{record.id}\n{protein}\n")

for r in SeqIO.parse("../data/proteins/spike_proteins.fasta", "fasta"):
    print(r.id, len(r.seq))
