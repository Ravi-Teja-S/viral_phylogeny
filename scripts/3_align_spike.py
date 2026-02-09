from pathlib import Path
import subprocess, shutil

BASE = Path(__file__).resolve().parent.parent
outdir = BASE / "data/aligned"
outdir.mkdir(exist_ok=True)

files = {
    BASE / "data/raw/spike_sequences.fasta": "spike_sequences.aln",
    BASE / "data/proteins/spike_proteins.fasta": "spike_proteins.aln"
}

for inp, aln_name in files.items():
    aln = outdir / aln_name
    dnd_raw = inp.with_suffix(".dnd")

    subprocess.run(
        ["clustalw", f"-INFILE={inp}", f"-OUTFILE={aln}", "-OUTPUT=CLUSTAL"],
        check=True
    )

    if dnd_raw.exists():
        shutil.move(dnd_raw, outdir / dnd_raw.name)

print("All alignments saved in:", outdir)
