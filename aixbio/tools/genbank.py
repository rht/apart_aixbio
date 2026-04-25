from __future__ import annotations

import logging
from io import StringIO

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

# pET-28a(+) backbone: 5369 bp.
# NOTE: The actual backbone sequence is not embedded due to licensing.
# The backbone is represented as N's with annotated functional elements
# at their approximate positions within the vector.
PET28A_BACKBONE_SIZE = 5369

# Approximate feature positions within pET-28a(+) (from SnapGene/Addgene maps)
_PET28A_FEATURES = [
    ("regulatory", 370, 386, {"note": ["T7 promoter"], "regulatory_class": ["promoter"]}),
    ("regulatory", 386, 404, {"note": ["lac operator"]}),
    ("regulatory", 404, 411, {"note": ["ribosome binding site"]}),
    ("CDS", 2970, 3785, {"gene": ["kan"], "note": ["Kanamycin resistance (KanR)"]}),
    ("rep_origin", 4338, 4927, {"note": ["pBR322 origin of replication"]}),
    ("rep_origin", 4990, 5446, {"note": ["f1 origin of replication"]}),
    ("terminator", 450, 497, {"note": ["T7 terminator"]}),
]


def build_plasmid_record(
    chain_id: str,
    cassette_dna: str,
    vector: str = "pET-28a(+)",
    cloning_sites: tuple[str, ...] = ("BamHI", "XhoI"),
) -> tuple[str, int]:
    insert_size = len(cassette_dna)
    total_size = PET28A_BACKBONE_SIZE + insert_size

    logger.warning(
        "GenBank backbone uses placeholder N's (not the real pET-28a(+) sequence). "
        "The output file is structurally valid but cannot be used for synthesis "
        "without replacing the backbone with the actual vector sequence."
    )

    plasmid_seq = "N" * PET28A_BACKBONE_SIZE + cassette_dna
    record = SeqRecord(
        Seq(plasmid_seq),
        id=f"{vector}_{chain_id}",
        name=f"{vector}_{chain_id}",
        description=f"{vector} with {chain_id} insert (WARNING: backbone is placeholder N's)",
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )

    # Annotate known functional elements on the backbone
    for feat_type, start, end, qualifiers in _PET28A_FEATURES:
        record.features.append(SeqFeature(
            FeatureLocation(start, end),
            type=feat_type,
            qualifiers=qualifiers,
        ))

    # Backbone source region
    record.features.append(SeqFeature(
        FeatureLocation(0, PET28A_BACKBONE_SIZE),
        type="source",
        qualifiers={"note": [f"{vector} backbone (placeholder — replace with actual sequence)"]},
    ))

    # Insert CDS
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
