from .BedmethylItem import BedmethylItem


class BedmethylItemGroup:
    """A group of BedmethylItem instances for reporting."""

    def __init__(self, bedmethylitems, name="Unnamed", methylated_cutoff=30):
        self.bedmethylitems = bedmethylitems
        self.name = name
        self.methylated_cutoff = methylated_cutoff  # see BedmethylItem class

    def perform_all_analysis_in_bedmethylitemgroup(self):
        for bedmethylitem in self.bedmethylitems:
            bedmethylitem.perform_analysis()
        self.comparisons_performed = True
