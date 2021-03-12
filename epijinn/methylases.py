from .epijinn import Methylase

EcoKDam = Methylase(name="EcoKDam", sequence="GATC", index_pos=1, index_neg=2)
EcoKDcm = Methylase(name="EcoKDcm", sequence="CCWGG", index_pos=1, index_neg=3)
EcoBI = Methylase(name="EcoBI", sequence="TGANNNNNNNNTGCT", index_pos=2, index_neg=11)
EcoKI = Methylase(name="EcoKI", sequence="AACNNNNNNGTGC", index_pos=1, index_neg=10)


METHYLASES = [
    EcoKDam,
    EcoKDcm,
    EcoBI,
    EcoKI,
]
