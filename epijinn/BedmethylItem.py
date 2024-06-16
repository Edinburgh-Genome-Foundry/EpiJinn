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

    def subset_bed_to_mod_subtype(self, mod):
        # Columnname from Bedmethyl specification
        bed_subtype = self.bed.loc[self.bed["modified_base_code_and_motif"] == mod]
        return bed_subtype

    def subset_bed_to_base_matches(self, modified_base="C"):
        mod_base = modified_base.upper()
        mod_base_complement = COMPLEMENTS[mod_base]
        # POSITIVE STRAND
        matching_positions_on_positive_strand = [
            pos for pos, base in enumerate(self.record) if base.upper() == mod_base
        ]  # for positive strand only
        positive_strand_filter = self.bed["start_position"].isin(
            matching_positions_on_positive_strand
        ) & (self.bed["strand"] == "+")

        # NEGATIVE STRAND
        matching_positions_on_negative_strand = [
            pos
            for pos, base in enumerate(self.record)
            if base.upper() == mod_base_complement
        ]
        negative_strand_filter = self.bed["start_position"].isin(
            matching_positions_on_negative_strand
        ) & (self.bed["strand"] == "-")

        bed_basematch = self.bed.loc[positive_strand_filter | negative_strand_filter]

        return bed_basematch
