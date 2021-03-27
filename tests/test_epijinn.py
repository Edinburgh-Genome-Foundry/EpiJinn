import epijinn

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def test_Methylase():
    EcoKDam = epijinn.Methylase(
        name="EcoKDam", sequence="GATC", index_pos=1, index_neg=2
    )
    assert EcoKDam.complement_table["R"] == "Y"
    assert EcoKDam.reverse("ACGT") == "TGCA"
    assert EcoKDam.complement("ACGT") == "TGCA"
    assert EcoKDam.reverse_complement("ACGT") == "ACGT"


def test_annotate_methylation():
    dna = Seq("TGACCCCCCCCTGCTCCCCCAGCACCCCCCCCTCA")
    dna_record = SeqRecord(dna, id="example", annotations={"molecule_type": "dna"})
    dna_annotated = epijinn.annotate_methylation(dna_record)
    assert len(dna_annotated.features) == 6
