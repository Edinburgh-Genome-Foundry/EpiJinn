# Copyright 2024 Edinburgh Genome Foundry, University of Edinburgh
#
# This file is part of EpiJinn.
#
# EpiJinn is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# EpiJinn is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Ediacara. If not, see <https://www.gnu.org/licenses/>.

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

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


class Methylator:
    """Class for finding methylation sites within a pattern (site) in a sequence.


    **Parameters**

    **sequence**
    > Sequence of ATGC characters (`str`).

    **site**
    > Sequence of restriction enzyme recognition site (`str`).

    **methylases**
    > Methylase class instances (`list`). Default uses built-in methylases.
    """

    def __init__(self, sequence, site, methylases=None):
        self.sequence = sequence
        if methylases is None:
            self.methylases = list(METHYLASES.values())
        else:
            self.methylases = methylases
        self.site = site

        self.pattern = dnachisel.SequencePattern(site)
        self.regions_seq = self.pattern.find_matches(self.sequence)

        self.site_rc = Methylase.reverse_complement(site)
        self.pattern_rc = dnachisel.SequencePattern(self.site_rc)
        self.regions_rc = self.pattern_rc.find_matches(self.sequence)

        self.regions = self.regions_seq + self.regions_rc

        self.report = ""

    def find_methylation_sites_in_pattern(self):
        """Run find_one_methylation_site_in_pattern() for each enzyme in methylases."""

        self.report += "Matches against methylase enzyme sites:\n\n"
        for methylase in self.methylases:
            self.find_one_methylation_site_in_pattern(methylase)
            self.report += "\n"

    def find_one_methylation_site_in_pattern(self, methylase):
        """Find overlapping methylation and restriction sites."""

        extended_regions = self.extend_restriction_regions(methylase)

        # For matching against positive strand of methylation pattern:
        expression = dnachisel.DnaNotationPattern.dna_sequence_to_regexpr(
            methylase.sequence
        )
        pattern = dnachisel.SequencePattern(expression)

        # For matching against negative strand of methylation pattern:
        expression_rc = dnachisel.DnaNotationPattern.dna_sequence_to_regexpr(
            methylase.rc
        )
        pattern_rc = dnachisel.SequencePattern(expression_rc)

        self.report += methylase.name + "\n"
        self.report += "=" * len(methylase.name) + "\n"

        for region in extended_regions:
            region_sequence = self.sequence[region.start : region.end]
            self.report += "Region:" + str(region) + "\n"

            match_location = pattern.find_matches(region_sequence)
            if len(match_location) != 0:
                self.report += "Match in positive strand: %s\n" % region_sequence
            else:
                self.report += "Positive strand: -\n"

            match_location_rc = pattern_rc.find_matches(region_sequence)
            if len(match_location_rc) != 0:
                self.report += "Match in negative strand: %s\n" % region_sequence
            else:
                self.report += "Negative strand: -\n"
            self.report += "\n"

    def extend_restriction_regions(self, methylase):
        """Modify list of dnachisel.Location of restriction sites to include
        flanking nucleotides around restriction sites.
        """

        extended_regions = []
        for region in self.regions:
            region = region.extended(
                methylase.index_neg, left=True, right=False
            )  # extension upstream

            m = len(methylase.sequence) - (methylase.index_pos + 1)
            region = region.extended(
                m, upper_limit=len(self.sequence), left=False, right=True
            )  # extension downstream
            extended_regions.append(region)
        return extended_regions


def annotate_methylation(seqrecord, methylases=None):
    """Annotate SeqRecord with methylation patterns and methylated nucleotides.

    **Parameters**

    **seqrecord**
    > A Biopython `SeqRecord`.

    **methylases**
    > A `list` of `Methylase` instances.
    """
    if methylases is None:
        methylases = list(METHYLASES.values())
    for methylase in methylases:
        name = methylase.name
        regex = dnachisel.DnaNotationPattern.dna_sequence_to_regexpr(methylase.sequence)
        pattern = dnachisel.SequencePattern(regex)
        match_location = pattern.find_matches(str(seqrecord.seq))
        if len(match_location) != 0:
            for match in match_location:
                label = "@epijinn(" + methylase.name + ")"
                seqrecord.features.append(
                    SeqFeature(
                        FeatureLocation(match.start, match.end),  # dnachisel.Location
                        type="misc_feature",
                        id="@epijinn",
                        qualifiers={"label": label, "note": name},
                    )
                )
                # Mark the methylation site for checking overlap with restriction site
                methylated_position = match.start + methylase.index_pos
                methylated_nucleotide = str(methylase.sequence[methylase.index_pos])
                label = "@epijinn(" + methylated_nucleotide + ", strand=1)"
                seqrecord.features.append(
                    SeqFeature(
                        FeatureLocation(
                            methylated_position, methylated_position + 1, strand=1
                        ),
                        type="misc_feature",
                        id="@epijinn",
                        qualifiers={"label": label, "note": name},
                    )
                )
                # Negative strand (don't annotate if pattern is 1-base long):
                if methylase.index_neg:  # `None` skips this step
                    methylated_position = match.start + methylase.index_neg
                    methylated_nucleotide = str(
                        Seq(
                            methylase.sequence[methylase.index_neg]
                        ).reverse_complement()
                    )
                    label = "@epijinn(" + methylated_nucleotide + ", strand=-1)"
                    seqrecord.features.append(
                        SeqFeature(
                            FeatureLocation(
                                methylated_position, methylated_position + 1, strand=-1
                            ),
                            type="misc_feature",
                            id="@epijinn",
                            qualifiers={"label": label, "note": name},
                        )
                    )

        # Repeat for reverse complement, if not palindromic:
        if Seq(methylase.sequence) != Seq(methylase.sequence).reverse_complement():
            name = methylase.name
            methylase_rc_seq = dnachisel.reverse_complement(methylase.sequence)
            regex = dnachisel.DnaNotationPattern.dna_sequence_to_regexpr(
                methylase_rc_seq
            )
            pattern = dnachisel.SequencePattern(regex)
            match_location = pattern.find_matches(str(seqrecord.seq))
            if len(match_location) != 0:
                for match in match_location:
                    label = "@epijinn_rc(" + methylase.name + ")"
                    seqrecord.features.append(
                        SeqFeature(
                            FeatureLocation(
                                match.start, match.end
                            ),  # dnachisel.Location
                            type="misc_feature",
                            id="@epijinn_rc",
                            qualifiers={"label": label, "note": name},
                        )
                    )
                    # Mark the methylation site:
                    # reverse complement so need to count backwards, and strand=-1
                    # subtract 1 to account for range
                    methylated_position = match.end - 1 - methylase.index_pos
                    methylated_nucleotide = str(methylase.sequence[methylase.index_pos])
                    label = "@epijinn(" + methylated_nucleotide + ", strand=-1)"
                    seqrecord.features.append(
                        SeqFeature(
                            FeatureLocation(
                                methylated_position, methylated_position + 1, strand=-1
                            ),
                            type="misc_feature",
                            id="@epijinn",
                            qualifiers={"label": label, "note": name},
                        )
                    )
                    # Mark methylation in antisense of the enzyme site:
                    # reverse complement so need to count backwards, and strand=1
                    # subtract 1 to account for range
                    if methylase.index_neg:  # again, don't annotate if 1-base long
                        methylated_position = match.end - 1 - methylase.index_neg
                        methylated_nucleotide = str(
                            Seq(
                                methylase.sequence[methylase.index_neg]
                            ).reverse_complement()
                        )
                        label = "@epijinn(" + methylated_nucleotide + ", strand=1)"
                        seqrecord.features.append(
                            SeqFeature(
                                FeatureLocation(
                                    methylated_position,
                                    methylated_position + 1,
                                    strand=1,
                                ),
                                type="misc_feature",
                                id="@epijinn",
                                qualifiers={"label": label, "note": name},
                            )
                        )

    return seqrecord


# Cytosine:
AluI = Methylase(name="AluI", sequence="AGCT", index_pos=2, index_neg=1)
BamHI = Methylase(name="BamHI", sequence="GGATCC", index_pos=4, index_neg=1)
CpG = Methylase(name="CpG", sequence="CG", index_pos=0, index_neg=1)
EcoKDcm = Methylase(name="EcoKDcm", sequence="CCWGG", index_pos=1, index_neg=3)
GpC = Methylase(name="GpC", sequence="GC", index_pos=1, index_neg=0)
HaeIII = Methylase(name="HaeIII", sequence="GGCC", index_pos=2, index_neg=1)
Hhal = Methylase(name="Hhal", sequence="GCGC", index_pos=1, index_neg=2)
HpaII = Methylase(name="HpaII", sequence="CCGG", index_pos=1, index_neg=2)
# For cytosine (not an enzyme):
MetC = Methylase(name="MetC", sequence="C", index_pos=0, index_neg=None)
MspI = Methylase(name="MspI", sequence="CCGG", index_pos=0, index_neg=3)

# Adenine:
EcoBI = Methylase(name="EcoBI", sequence="TGANNNNNNNNTGCT", index_pos=2, index_neg=11)
EcoGII = Methylase(name="EcoGII", sequence="A", index_pos=0, index_neg=None)
EcoKDam = Methylase(name="EcoKDam", sequence="GATC", index_pos=1, index_neg=2)
EcoKI = Methylase(name="EcoKI", sequence="AACNNNNNNGTGC", index_pos=1, index_neg=10)
EcoRI = Methylase(name="EcoRI", sequence="GAATTC", index_pos=2, index_neg=3)
TaqI = Methylase(name="TaqI", sequence="TCGA", index_pos=3, index_neg=0)


METHYLASES = {
    AluI.name: AluI,
    BamHI.name: BamHI,
    CpG.name: CpG,
    EcoKDcm.name: EcoKDcm,
    GpC.name: GpC,
    HaeIII.name: HaeIII,
    Hhal.name: Hhal,
    HpaII.name: HpaII,
    MetC.name: MetC,
    MspI.name: MspI,
    EcoBI.name: EcoBI,
    # EcoGII.name: EcoGII,
    EcoKDam.name: EcoKDam,
    EcoKI.name: EcoKI,
    EcoRI.name: EcoRI,
    TaqI.name: TaqI,
}
