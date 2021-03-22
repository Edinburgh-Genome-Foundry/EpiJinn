<p align="center">
<img alt="EpiJinn logo" title="EpiJinn" src="images/epijinn.png" width="120">
</p>


# EpiJinn

[![Build Status](https://travis-ci.org/Edinburgh-Genome-Foundry/EpiJinn.svg?branch=main)](https://travis-ci.org/Edinburgh-Genome-Foundry/EpiJinn)
[![Coverage Status](https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/EpiJinn/badge.svg?branch=main)](https://coveralls.io/github/Edinburgh-Genome-Foundry/EpiJinn?branch=main)

**Work in progress!**

**EpiJinn** checks whether recognition sites of prokaryotic DNA *methylase enzymes* overlap with a recognition site of a *restriction enzyme*, in a DNA sequence. Methylation of restriction site nucleotides blocks recognition/restriction and thus DNA assembly.

The module contains the `Methylator` class for storing a sequence, methylation enzymes and a restriction enzyme recognition site. It has a method for finding overlaps, and uses [DNA Chisel](https://edinburgh-genome-foundry.github.io/DnaChisel/) to find sequence matches.

An example overlap:

    ...ccgcatgaagggcgcgccaggtctcaccctgaattcgcg...
                          ggtctc    : BsaI restriction site
                       CCAGGTCTCACC : Match in positive strand
                       CCWGG        : Dcm methylation site
                        *           : methylated cytosine
                          *         : methylated cytosine (on other strand)

For information on the effect of DNA methylation on each enzyme, see the [Restriction Enzyme Database](http://rebase.neb.com/rebase/rebms.html).


## Usage

```python
import epijinn
methylator = epijinn.Methylator(sequence=str(sequence), site=site_BsaI)
methyl.find_methylation_sites_in_pattern()
```


## Example

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


## DNA sulfur modification

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
