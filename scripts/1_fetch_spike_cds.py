import time
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent

Entrez.email = "your_email@example.com"
DELAY = 1.2

ACCESSIONS = {
    "Wuhan": "NC_045512.2",

    "Alpha_B117": "MW422255",
    "Beta_B1351": "MW598419",
    "Gamma_P1": "MT126808",
    "Delta_B16172": "MZ208926",

    "Omicron_BA1": "OM570283",
    "Omicron_BA2": "OM617939",
    "Omicron_BA5": "OP093341",

    "Lambda_C37": "MZ068159",
    "Mu_B1621": "MZ344997",

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


from Bio.SeqFeature import CompoundLocation

def extract_spike(record):
    for feature in record.features:
        if feature.type != "CDS":
            continue

        text = " ".join(
            feature.qualifiers.get("gene", []) +
            feature.qualifiers.get("product", [])
        ).lower()

        if "spike" in text or "surface glycoprotein" in text:
            loc = feature.location

            # Handle joined/multi-part CDS correctly
            if isinstance(loc, CompoundLocation):
                seq = ""
                for part in loc.parts:
                    seq += part.extract(record.seq)
                return seq

            return loc.extract(record.seq)

    raise ValueError(f"Spike CDS not found in {record.id}")


with open(OUTPUT, "w") as out:
    for name, acc in ACCESSIONS.items():
        print("Fetching", name)

        genome = fetch_genbank(acc)
        spike = extract_spike(genome)

        print(f"{name:15} spike length (nt): {len(spike)}")

        rec = SeqRecord(
            spike,
            id=name,
            description=f"Spike CDS {acc}"
        )

        SeqIO.write(rec, out, "fasta")
        time.sleep(DELAY)
