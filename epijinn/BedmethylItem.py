import matplotlib.pyplot as plt

import dna_features_viewer


class BedmethylItem:
    """Analyse a bedmethyl file.


    **Parameters**

    **reference**
    > A Biopython SeqRecord.

    **bedmethyl**
    > A pandas dataframe of a bedmethyl file for the reference.
    """

    def __init__(self, reference, bedmethyl):
        self.record = reference
        self.name = reference.name
        self.id = reference.id
        self.bed = bedmethyl
        self.methylated_cutoff = 0.6  # fraction of methylated base at a position
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
