import time
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent

Entrez.email = "your_email@example.com"
DELAY = 1.2

ACCESSIONS = {
    "Wuhan": "NC_045512.2",
    "Delta": "MZ208926",
    "Omicron": "OM570283",
    "Bat_outgroup": "MN996532"
}

target_dir = BASE / "data/raw"
target_dir.mkdir(parents=True,exist_ok=True)

OUTPUT = target_dir/ "spike_sequences.fasta"


def fetch_genbank(acc):
    handle = Entrez.efetch(
        db="nucleotide",
        id=acc,
        rettype="gb",
        retmode="text"
    )
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record


def extract_spike(record):
    for feature in record.features:
        if feature.type == "CDS":
            product = feature.qualifiers.get("product", [""])[0].lower()
            gene = feature.qualifiers.get("gene", [""])[0].lower()
            if "spike" in product or gene == "s":
                return feature.extract(record.seq)
    return None


with open(OUTPUT, "w") as out:
    for name, acc in ACCESSIONS.items():
        print("Fetching", name)
        genome = fetch_genbank(acc)
        spike = extract_spike(genome)

        rec = SeqRecord(
            spike,
            id=name,
            description=f"Spike CDS {acc}"
        )

        SeqIO.write(rec, out, "fasta")
        time.sleep(DELAY)

print("Saved →", OUTPUT)
