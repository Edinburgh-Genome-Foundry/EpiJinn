# Copyright 2024 Edinburgh Genome Foundry, University of Edinburgh
#
# This file is part of EpiJinn.
#
# EpiJinn is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# EpiJinn is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with Ediacara. If not, see <https://www.gnu.org/licenses/>.

import os

import pandas

import Bio

from .BedmethylItem import BedmethylItem

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

    def perform_all_analysis_in_bedmethylitemgroup(self):
        for bedmethylitem in self.bedmethylitems:
            bedmethylitem.perform_analysis()
        self.comparisons_performed = True

    @staticmethod
    def subset_bed_columns(bed):
        bed_report = bed[columns_for_pdf_report]
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
