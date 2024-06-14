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

MODIFICATIONS = {
    "h": "C",
    "m": "C",
}


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
