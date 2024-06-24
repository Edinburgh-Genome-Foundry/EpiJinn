<p align="center">
<img alt="EpiJinn logo" title="EpiJinn" src="images/epijinn.png" width="120">
</p>

# EpiJinn

[![Build Status](https://travis-ci.org/Edinburgh-Genome-Foundry/EpiJinn.svg?branch=main)](https://travis-ci.org/Edinburgh-Genome-Foundry/EpiJinn)
[![Coverage Status](https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/EpiJinn/badge.svg?branch=main)](https://coveralls.io/github/Edinburgh-Genome-Foundry/EpiJinn?branch=main)

**Work in progress!**

**EpiJinn** is a Python package for working with modified (methylated) nucleotides.

* Create a readable report from *bedMethyl* files created by [modkit](https://github.com/nanoporetech/modkit/).
* Annotate prokaryotic DNA *methylase enzyme* recognition sites in a [Biopython](https://biopython.org/) SeqRecord.
* Check whether recognition sites of prokaryotic DNA methylase enzymes *overlap* with a recognition site of a *restriction enzyme*, in a DNA sequence. Methylation of restriction site nucleotides blocks recognition/restriction and thus DNA assembly.

The software is geared towards working with plasmids. Several more future functionalities are *planned:* comparing methylation status with expected methylation levels; methylase site recognition etc.

## Install

```bash
pip install git+https://github.com/Edinburgh-Genome-Foundry/EpiJinn.git
```

See additional install instructions for the [PDF Reports](https://github.com/Edinburgh-Genome-Foundry/pdf_reports) dependency,
and its Weasyprint dependency.

## Usage

### bedMethyl files

```python
import epijinn
bedmethylitemgroup = epijinn.read_sample_sheet(
    sample_sheet="sample_sheet.csv",
    genbank_dir='genbank',
    bedmethyl_dir='bedmethyl',
    parameter_sheet='param_sheet.csv',)
bedmethylitemgroup.perform_all_analysis_in_bedmethylitemgroup()
epijinn.write_bedmethylitemgroup_report(bedmethylitemgroup=bedmethylitemgroup, pdf_file="report.pdf", html_file="report.html")
```

Both `pdf_file` and `html_file` are optional, specify `None` to exclude either of them.
An example sample sheet and parameter sheet is included in the `examples` directory.

### Recognition site overlap

The module contains the `Methylator` class for storing a sequence, methylation enzymes and a restriction enzyme recognition site. It has a method for finding overlaps, and uses [DNA Chisel](https://edinburgh-genome-foundry.github.io/DnaChisel/) to find sequence matches.

An example overlap:

    ...ccgcatgaagggcgcgccaggtctcaccctgaattcgcg...
                          ggtctc    : BsaI restriction site
                       CCAGGTCTCACC : Match in positive strand
                       CCWGG        : Dcm methylation site
                        *           : methylated cytosine
                          *         : methylated cytosine (on other strand)

For information on the effect of DNA methylation on each enzyme, see the [Restriction Enzyme Database](http://rebase.neb.com/rebase/rebms.html).

```python
import epijinn
methylator = epijinn.Methylator(sequence=str(sequence), site=site_BsaI)
methyl.find_methylation_sites_in_pattern()
```

### Example

```python
import epijinn
import Bio

sequence = 'ATGTCCCCATGCCTAC' + 'AGCAAGGC' + 'CGTCTC' + 'A' + 'GGCCCCCCCCCCCCA'  # seq + EcoBI (+ BsmBI +) EcoBI + seq

rest_dict = Bio.Restriction.Restriction_Dictionary.rest_dict
site_BsmBI = rest_dict['BsmBI']['site']

epijinn.EcoBI.sequence
# 'TGANNNNNNNNTGCT'
methylator = epijinn.Methylator(sequence, site=site_BsmBI)
methylator.find_methylation_sites_in_pattern()
print(methylator.report)
```

Result:

    Matches against methylase enzyme sites:

    EcoKDam
    =======
    Region: 22-32(+)
    Positive strand: -
    Negative strand: -


    EcoKDcm
    =======
    Region: 21-33(+)
    Positive strand: -
    Negative strand: -


    EcoBI
    =====
    Region: 13-42(+)
    Positive strand: -
    Match in negative strand: TACAGCAAATCCGTCTCAGGCCCCCCCCC


    EcoKI
    =====
    Region: 14-41(+)
    Positive strand: -
    Negative strand: -

### DNA sulfur modification

The same approach can be used for finding enzyme site overlaps with other epigenetic modifications. For example, in DNA phosphorothioation, an oxygen on the DNA backbone is replaced with sulfur.

```python
thio = epijinn.Methylator(sequence, site=site_BsmBI, methylases=epijinn.DND)
thio.find_methylation_sites_in_pattern()
```

This returns an overlap with a putative *dnd* target site of *Streptomyces lividans 1326* with conserved sequence GGCC:

    Dnd_Sli1326
    ===========
    Region: 21-33(+)
    Match in positive strand: GGCCGTCTCAGG
    Match in negative strand: GGCCGTCTCAGG

## Versioning

EpiJinn uses the [semantic versioning](https://semver.org) scheme.

## License = GPLv3+

Copyright 2024 Edinburgh Genome Foundry, University of Edinburgh

EpiJinn was written at the [Edinburgh Genome Foundry](https://edinburgh-genome-foundry.github.io/)
by [Peter Vegh](https://github.com/veghp), and is released under the GPLv3 license.
