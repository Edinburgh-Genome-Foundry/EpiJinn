# Copyright 2024 Edinburgh Genome Foundry, University of Edinburgh
#
# This file is part of EpiJinn.
#
# EpiJinn is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# EpiJinn is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Ediacara. If not, see <https://www.gnu.org/licenses/>.

import copy
import os

import pandas
import matplotlib.pyplot as plt

import Bio

import dna_features_viewer

from .Methyl import METHYLASES
from .Methyl import annotate_methylation

COMPLEMENTS = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
}

# From https://github.com/nanoporetech/modkit/blob/master/book/src/intro_bedmethyl.md
BEDMETHYL_HEADER = [
    "chrom",
    "start_position",
    "end_position",
    "modified_base_code_and_motif",
    "score",
    "strand",
    "strand_start_position",
    "strand_end_position",
    "color",
    "Nvalid_cov",
    "percent_modified",
    "Nmod",
    "Ncanonical",
    "Nother_mod",
    "Ndelete",
    "Nfail",
    "Ndiff",
    "Nnocall",
]

# Remove duplicate and unnecessary columns from report:
columns_for_pdf_report = [
    BEDMETHYL_HEADER[i] for i in [1, 5, 9, 10, 11, 12, 13, 14, 15, 16, 17]
] + [
    "status"
]  # added during binarization

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")
# For looking up modified_base_code_and_motif entries in bedmethyl:
MODIFICATION_CODES = pandas.read_csv(os.path.join(DATA_DIR, "mod_base_codes.csv"))
# Adapted From https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf


class CustomTranslator(dna_features_viewer.BiopythonTranslator):
    """Custom translator used for the main plot."""

    def compute_filtered_features(self, features):
        """Display only selected features."""
        filtered_features = []
        n = 8  # a good number of features to display
        # Keep longest features
        if len(features) > n:
            feature_lengths = []
            for feature in features:
                feature_lengths += [len(feature.location)]
            feature_lengths.sort(reverse=True)
            max_length = feature_lengths[n]
            for feature in features:
                if len(feature.location) > max_length:
                    filtered_features += [feature]
        else:  # no need to do anything if not enough features
            filtered_features = features

        return filtered_features


class PatternTranslator(dna_features_viewer.BiopythonTranslator):
    """Custom translator used for the pattern annotation plot."""

    def compute_feature_color(self, feature):
        status = "status"
        # Values are from binarize_bed() :
        if feature.qualifiers[status] == "1":
            return "red"
        elif feature.qualifiers[status] == "U":
            return "yellow"
        elif feature.qualifiers[status] == "0":
            return "grey"
        # to be implemented later:
        elif feature.qualifiers[status] == "?":
            return "tab:cyan"
        # There should be no other options.


class BedResult:
    """Results of a bedmethyl table analysis."""

    def __init__(self, modification, methylase, bed, record):
        self.feature_cutoff = 50  # do not create a plot if there are more features
        self.modification = modification
        self.mod_abbreviation = MODIFICATION_CODES.loc[
            MODIFICATION_CODES["Code"] == self.modification
        ].Abbreviation.iloc[
            0
        ]  # there is only one entry so we take the first
        self.mod_name = MODIFICATION_CODES.loc[
            MODIFICATION_CODES["Code"] == self.modification
        ].Name.iloc[0]
        self.mod_chebi = MODIFICATION_CODES.loc[
            MODIFICATION_CODES["Code"] == self.modification
        ].ChEBI.iloc[0]
        self.unmodified_base = MODIFICATION_CODES.loc[
            MODIFICATION_CODES["Code"] == self.modification
        ].Unmodified_base.iloc[0]

        self.methylase = methylase
        self.methylase_str = methylase.name  # for the report
        self.methylase_seq = methylase.sequence  # for the report
        self.bed = bed
        self.record = record

        # FILTER (ANNOTATED) RECORD AND ADD STATUS:
        filtered_features = []
        labelstart = (
            "@epijinn(" + self.unmodified_base + ", strand="
        )  # unfinished to account for +/- strands
        for feature in record.features:
            if feature.id == "@epijinn":  # as annotated by the function
                if feature.qualifiers["label"].startswith(labelstart):
                    # To display seq in report:
                    feature.qualifiers["label"] = self.unmodified_base
                    # Assign status:
                    # no need to check for strand as the complement won't be modified
                    location_start = int(feature.location.start)
                    # expect exactly one bed table entry per modified nucleotide:
                    status_symbol = self.bed.loc[
                        self.bed["LOC"] == location_start
                    ].STATUS.iloc[0]
                    # used by a custom BiopythonTranslator to colour the annotation:
                    feature.qualifiers["status"] = status_symbol
                    filtered_features += [feature]
        record.features = filtered_features

        if len(self.record.features) > self.feature_cutoff:
            self.img_created = False
        else:
            self.img_created = True
            fig, ax1 = plt.subplots(1, 1, figsize=(8, 3))
            graphic_record = PatternTranslator().translate_record(self.record)
            graphic_record.plot(ax=ax1, with_ruler=False, strand_in_label_threshold=4)

            self.plot = fig


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

    def perform_analysis(self, parameter_dict):
        """Perform analysis and plot the sequence."""
        self.fig = self.plot_record()
        self.results = []  # BedResult instances. For easy reference in pug template.
        methylase_strings = parameter_dict["methylases"].split(" ")
        # space-separated as per specifications
        for methylase_str in methylase_strings:
            methylase = METHYLASES[methylase_str]
            for modification in self.bed.modified_base_code_and_motif.unique():
                # RUN BED SUBSET ETC
                bed_subtype = self.subset_bed_to_mod_subtype(self.bed, mod=modification)
                annotated_record, bed_pattern_match = self.subset_bed_to_pattern_match(
                    methylase, bed=bed_subtype
                )
                bed_binarized = self.binarize_bed(
                    bed_pattern_match,
                    met_cutoff=parameter_dict["methylated_cutoff"],
                    nonmet_cutoff=parameter_dict["unmethylated_cutoff"],
                )
                final_bed = self.subset_bed_columns(bed_binarized)
                bedresult = BedResult(
                    modification=modification,
                    methylase=methylase,
                    bed=final_bed,
                    record=annotated_record,
                )
                self.results += [bedresult]

    def annotate_record(self):
        return self.record

    def plot_record(self):
        fig, ax1 = plt.subplots(1, 1, figsize=(8, 2))
        graphic_record = CustomTranslator().translate_record(self.record)
        graphic_record.plot(ax=ax1, with_ruler=True, strand_in_label_threshold=4)

        return fig

    @staticmethod
    def subset_bed_columns(bed):
        bed_report = bed[columns_for_pdf_report]
        # These were designed to be more informative and fit the report:
        new_columnnames = [
            "LOC",
            "Strand",
            "COV",
            "% mod",
            "MOD",
            "STD",
            "OTH",
            "del",
            "fail",
            "diff",
            "nocall",
            "STATUS",
        ]
        bed_report.columns = new_columnnames

        return bed_report

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
        record_copy = copy.deepcopy(self.record)
        annotated_record = annotate_methylation(record_copy, methylases=[methylase])

        # CREATE LIST OF POSITIONS
        positive_strand_locations = []
        negative_strand_locations = []
        label_pos = "@epijinn(" + mod_base + ", strand=1)"
        label_neg = "@epijinn(" + mod_base + ", strand=-1)"
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
    def binarize_bed(bed, met_cutoff, nonmet_cutoff):
        bed["status"] = "U"  # prefill undetermined
        for index, row in bed.iterrows():
            if row["percent_modified"] >= float(met_cutoff) * 100:
                bed.loc[index, "status"] = "1"  # methylated, symbol may change
            elif row["percent_modified"] <= float(nonmet_cutoff) * 100:
                bed.loc[index, "status"] = "0"  # unmethylated
            # symbol "?" reserved for low coverage to be implemented later
        return bed


class BedmethylItemGroup:
    """A group of BedmethylItem instances for reporting.


    **Parameters**

    **bedmethylitems**
    > A list of BedmethylItem instances.

    **parameter_dict**
    > A dictionary of analysis parameters (`dict`).
    """

    def __init__(self, bedmethylitems, parameter_dict):
        self.bedmethylitems = bedmethylitems
        self.parameter_dict = parameter_dict
        self.number_of_samples = len(bedmethylitems)

    def perform_all_analysis_in_bedmethylitemgroup(self):
        for bedmethylitem in self.bedmethylitems:
            bedmethylitem.perform_analysis(parameter_dict=self.parameter_dict)
        self.comparisons_performed = True


def read_sample_sheet(
    sample_sheet, genbank_dir="", bedmethyl_dir="", parameter_sheet=""
):
    """Read a sample sheet into a BedmethylItemGroup.


    **Parameters**

    **sample_sheet**
    > CSV file path (`str`). No header and columns must be in this order: projectname,
    sample, Genbank name (without extension), bedmethyl file.

    **genbank_dir**
    > Directory of the Genbank files (`str`). Default: local directory.

    **bedmethyl_dir**
    > Directory of the bedmethyl files (`str`). Default: local directory.

    **parameter_sheet**
    > CSV file path (`str`). Use 'Parameter', 'Value' header for columns. If a
    'projectname' is specified, it overwrites the sample sheet value.
    """
    # READ PARAMETERS
    param_df = pandas.read_csv(parameter_sheet, usecols=["Parameter", "Value"])
    parameter_dict = dict(param_df.values)

    # READ SAMPLES
    sample_df = pandas.read_csv(sample_sheet, header=None)
    # add columnnames

    # CREATE ITEMS
    # We allow Sequeduct to specify the projectname as a command parameter as well;
    if not "projectname" in parameter_dict:
        # first entry of the first column (contains projectname):
        parameter_dict["projectname"] = sample_df.iloc[0, 0]
    # Allow the user to not specify these in the sheet:
    # The default cutoffs are based on preliminary data.
    if not "methylated_cutoff" in parameter_dict:
        parameter_dict["methylated_cutoff"] = 0.7
    if not "unmethylated_cutoff" in parameter_dict:
        parameter_dict["unmethylated_cutoff"] = 0.3

    bedmethylitems = []
    for index, row in sample_df.iterrows():
        genbank_name = row[2]  # number specified by the sample sheet format
        genbank_path = os.path.join(genbank_dir, genbank_name + ".gb")  # Genbank ext
        record = Bio.SeqIO.read(genbank_path, "genbank")
        record.id = genbank_name
        record.name = genbank_name
        record.annotations["molecule_type"] = "DNA"

        bed_name = row[3]  # number specified by the sample sheet format
        bed_path = os.path.join(bedmethyl_dir, bed_name)
        bed_df = pandas.read_csv(bed_path, header=None, delimiter="\t")
        bed_df.columns = BEDMETHYL_HEADER
        bedmethylitems += [
            BedmethylItem(sample=row[1], reference=record, bedmethyl=bed_df)
        ]  # number specified by the sample sheet format

    bedmethylitemgroup = BedmethylItemGroup(
        bedmethylitems=bedmethylitems, parameter_dict=parameter_dict
    )
    return bedmethylitemgroup
