import epijinn


def test_Methylase():
    EcoKDam = epijinn.Methylase(
        name="EcoKDam", sequence="GATC", index_pos=1, index_neg=2
    )
    assert EcoKDam.reverse("ACGT") == "TGCA"
