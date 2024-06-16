# Copyright 2024 Edinburgh Genome Foundry, University of Edinburgh
#
# This file is part of EpiJinn.
#
# EpiJinn is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# EpiJinn is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Ediacara. If not, see <https://www.gnu.org/licenses/>.

import matplotlib.pyplot as plt

import dna_features_viewer

from .epijinn import annotate_methylation

COMPLEMENTS = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
}


class BedmethylItem:
    """Analyse a bedmethyl file.


    **Parameters**

    **sample**
    > The name of the sample, for example a barcode id (`str`).

    **reference**
    > A Biopython SeqRecord.

    **bedmethyl**
    > A pandas dataframe of a bedmethyl file for the reference.
    """

    def __init__(self, sample, reference, bedmethyl):
        self.sample = sample
        self.record = reference
        self.name = reference.name
        self.id = reference.id
        self.bed = bedmethyl
        self.reference_length = len(self.record)

    def perform_analysis(self):
        """Perform analysis and plot the sequence."""
        self.annotated_record = self.annotate_record()
        self.fig = self.plot_record()

    def annotate_record(self):
        return self.record

    def plot_record(self):
        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(7, 3), sharex=True, gridspec_kw={"height_ratios": [4, 1]}
        )
        graphic_record = dna_features_viewer.BiopythonTranslator().translate_record(
            self.record
        )
        graphic_record.plot(ax=ax1, with_ruler=False, strand_in_label_threshold=4)

        return fig

    @staticmethod
    def subset_bed_to_mod_subtype(bed, mod):
        # Columnname from Bedmethyl file specification
        bed_subtype = bed.loc[bed["modified_base_code_and_motif"] == mod]

        return bed_subtype

    def subset_bed_to_base_matches(self, bed=None, modified_base="C"):
        if bed is None:  # optional bed allows linking multiple bed subset methods
            bed = self.bed
        mod_base = modified_base.upper()
        mod_base_complement = COMPLEMENTS[mod_base]
        # POSITIVE STRAND
        matching_positions_on_positive_strand = [
            pos
            for pos, base in enumerate(str(self.record.seq.upper()))
            if base == mod_base
        ]
        # Columnname from Bedmethyl file specification:
        positive_strand_filter = bed["start_position"].isin(
            matching_positions_on_positive_strand
        ) & (bed["strand"] == "+")

        # NEGATIVE STRAND
        matching_positions_on_negative_strand = [
            pos
            for pos, base in enumerate(str(self.record.seq.upper()))
            if base == mod_base_complement
        ]
        negative_strand_filter = bed["start_position"].isin(
            matching_positions_on_negative_strand
        ) & (bed["strand"] == "-")

        bed_basematch = bed.loc[positive_strand_filter | negative_strand_filter]

        return bed_basematch

    def subset_bed_to_pattern_match(self, methylase, bed=None):
        if bed is None:  # optional bed allows linking multiple bed subset methods
            bed = self.bed
        methylated_index = methylase.index_pos
        mod_base = methylase.sequence[methylated_index].upper()

        annotated_record = annotate_methylation(self.record)

        # CREATE LIST OF POSITIONS
        positive_strand_locations = []
        negative_strand_locations = []
        label_pos = "@epijinn(met" + mod_base + ", strand=1)"
        label_neg = "@epijinn(met" + mod_base + ", strand=-1)"
        for feature in annotated_record.features:
            if feature.id == "@epijinn":  # as annotated by function above
                if feature.qualifiers["label"] == label_pos:
                    positive_strand_locations += [feature.location.start]
                elif feature.qualifiers["label"] == label_neg:
                    negative_strand_locations += [feature.location.start]

        # SUBSET USING POSITIONS
        # Columnname from Bedmethyl file specification:
        positive_strand_filter = bed["start_position"].isin(
            positive_strand_locations
        ) & (bed["strand"] == "+")

        negative_strand_filter = bed["start_position"].isin(
            negative_strand_locations
        ) & (bed["strand"] == "-")

        bed_pattern_match = bed.loc[positive_strand_filter | negative_strand_filter]

        return annotated_record, bed_pattern_match

    @staticmethod
    def binarize_bed(bed):
        bed["status"] = "U"  # prefill undetermined
        for index, row in bed.iterrows():
            if row["percent_modified"] >= 70:  # good cutoff based on literature
                bed.loc[index, "status"] = "1"  # methylated, symbol may change
            elif row["percent_modified"] <= 30:  # cutoff for non-methylated
                bed.loc[index, "status"] = "0"  # unmethylated
            # symbol "?" reserved for low coverage to be implemented later
        return bed
