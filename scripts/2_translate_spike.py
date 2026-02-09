from Bio import SeqIO

with open("../data/proteins/spike_proteins.fasta", "w") as out:
    for record in SeqIO.parse("../data/raw/spike_sequences.fasta", "fasta"):
        protein = record.seq.translate(to_stop=True)
        out.write(f">{record.id}\n{protein}\n")
