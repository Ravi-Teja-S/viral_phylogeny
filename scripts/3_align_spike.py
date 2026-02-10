from pathlib import Path
import subprocess

BASE = Path(__file__).resolve().parent.parent
outdir = BASE / "data/aligned"
outdir.mkdir(exist_ok=True)

files = {
    BASE / "data/raw/spike_sequences.fasta": "spike_sequences.aln",
    BASE / "data/proteins/spike_proteins.fasta": "spike_proteins.aln"
}

for inp, aln_name in files.items():
    aln = outdir / aln_name

    print("Aligning:", inp.name)

    subprocess.run(
        ["mafft", "--auto", str(inp)],
        stdout=open(aln, "w"),
        check=True
    )

print("\nAll MAFFT alignments saved in:", outdir)
