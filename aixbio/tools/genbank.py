from __future__ import annotations

from io import StringIO

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


PET28A_BACKBONE_SIZE = 5369


def build_plasmid_record(
    chain_id: str,
    cassette_dna: str,
    vector: str = "pET-28a(+)",
    cloning_sites: tuple[str, ...] = ("BamHI", "XhoI"),
) -> tuple[str, int]:
    insert_size = len(cassette_dna)
    total_size = PET28A_BACKBONE_SIZE + insert_size

    plasmid_seq = "N" * PET28A_BACKBONE_SIZE + cassette_dna
    record = SeqRecord(
        Seq(plasmid_seq),
        id=f"{vector}_{chain_id}",
        name=f"{vector}_{chain_id}",
        description=f"{vector} with {chain_id} insert",
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )

    record.features.append(SeqFeature(
        FeatureLocation(0, PET28A_BACKBONE_SIZE),
        type="source",
        qualifiers={"note": [f"{vector} backbone"]},
    ))
    record.features.append(SeqFeature(
        FeatureLocation(PET28A_BACKBONE_SIZE, total_size),
        type="CDS",
        qualifiers={
            "gene": [chain_id],
            "note": [f"Optimized {chain_id} expression cassette"],
            "codon_start": [1],
        },
    ))

    buf = StringIO()
    from Bio import SeqIO
    SeqIO.write(record, buf, "genbank")
    return buf.getvalue(), insert_size
