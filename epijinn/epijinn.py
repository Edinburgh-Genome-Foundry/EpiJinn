import Bio
import dnachisel


class Methylase:
    """Methylase enzyme class.


    **Parameters**

    **name**
    > Name of the enzyme (`str`).

    **sequence**
    > Sequence of extended nucleotide characters (`str`).

    **index_pos**
    > Index of the methylated base on the positive strand (`int`).

    **index_neg**
    > Index of the methylated base on the negative strand (`int`).
    """

    complement_table = {
        "A": "T",
        "G": "C",
        "C": "G",
        "T": "A",
        "Y": "R",
        "R": "Y",
        "W": "W",
        "S": "S",
        "K": "M",
        "M": "K",
        "D": "H",
        "V": "B",
        "H": "D",
        "B": "V",
        "X": "X",
        "N": "N",
    }

    def __init__(self, name, sequence, index_pos, index_neg):
        self.name = name
        self.sequence = sequence
        self.rc = self.reverse_complement(sequence)
        self.index_pos = index_pos
        self.index_neg = index_neg

    @staticmethod
    def reverse(sequence):
        reverse = sequence[::-1]
        return reverse

    @staticmethod
    def complement(sequence):
        complement_letters = [Methylase.complement_table[letter] for letter in sequence]
        complement = "".join(complement_letters)
        return complement

    @staticmethod
    def reverse_complement(sequence):
        r = Methylase.reverse(sequence)
        rc = Methylase.complement(r)
        return rc
