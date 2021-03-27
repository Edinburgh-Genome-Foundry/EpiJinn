import epijinn

import Bio
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


def test_Methylator():
    rest_dict = Bio.Restriction.Restriction_Dictionary.rest_dict

    site_BsmBI = rest_dict["BsmBI"]["site"]

    methylator = epijinn.Methylator(
        "ATGTCCCCATGCCTACAGCAAGGCCGTCTCAGGCCCCCCCCCCCCA", site=site_BsmBI
    )

    methylator.find_methylation_sites_in_pattern()
    assert methylator.report


def test_annotate_methylation():
    dna = Seq("TGACCCCCCCCTGCTCCCCCAGCACCCCCCCCTCA")
    dna_record = SeqRecord(dna, id="example", annotations={"molecule_type": "dna"})
    dna_annotated = epijinn.annotate_methylation(dna_record)
    assert len(dna_annotated.features) == 6
