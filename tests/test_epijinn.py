import epijinn


def test_Methylase():
    EcoKDam = epijinn.Methylase(
        name="EcoKDam", sequence="GATC", index_pos=1, index_neg=2
    )
    assert EcoKDam.complement_table["R"] == "Y"
    assert EcoKDam.reverse("ACGT") == "TGCA"
    assert EcoKDam.complement("ACGT") == "TGCA"
    assert EcoKDam.reverse_complement("ACGT") == "ACGT"
