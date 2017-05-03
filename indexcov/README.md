indexcov
========

Quickly estimate coverage from a *whole-genome* bam or [**cram**](#CRAM) index. 
A bam index has 16KB resolution so that's what this gives, but it provides what appears to be a high-quality 
coverage estimate in seconds per genome.

The output is scaled to around 1. So a long stretch with values of 1.5 would be a heterozygous duplication.
This is useful as a quick QC to get coverage values across the genome.

In our tests, we can **estimate depth across 60X genomes for 30 samples in 30 seconds**.

Interactive HTML plots of depth are output for each chromosome. **Live examples of the interactive output are available [here](http://indexcov.s3-website-us-east-1.amazonaws.com/)**

Usage
=====

```
goleft indexcov --directory my-project-dir/ *.bam
```

This will create a number of text files described in the [Files](#Files) section below.

In addition, it will write a few `.html` files containing interactive plots.

For example, if we view the $prefix-indexcov-depth-X.html file for **X chromosome** we can see a
nice separation of samples by sex except at the PAR at the left:

![X Example](https://cloud.githubusercontent.com/assets/1739/21597648/074f06ca-d10b-11e6-8732-e9a2e8d1ecb5.png "x example")

That plot is taken directly from the HTML output by `indexcov`.

Using that separation, `indexcov` infers the copy-number of the sex chromosomes, outputs a stub .ped/.fam file with that
information, and makes a plot like this one:

![Sex Example](https://cloud.githubusercontent.com/assets/1739/21627994/2973d464-d1d9-11e6-9962-5d3ac0f80329.png "sex example")
Where here the males and females separate by the X and Y chromosomes perfectly.

In some cases, we have found *XXY* and *XYY* samples this way.


`indexcov` will output a coverage (ROC) plot that shows how much of the genome is coverage at at given (scaled) depth.
This is output to a $prefix-depth-roc.html file and looks like:
![ROC Example](https://cloud.githubusercontent.com/assets/1739/21599983/b27fa4d8-d132-11e6-95b9-e9fa8ae64412.png "ROC example")

Here we can see that one sample has much lower coverage than the rest, and we can hover and determine the exact sample.


Finally, `indexcov` will output a `$prefix-indexcov-bins.html` file with a point per sample. Samples with high
values on the y-axis have very uneven coverage (this will affect SV calling). Samples with high values on There
x-axis have many missing bins (likely truncated bam files).

<a name="CRAM"></a> CRAM
========================

CRAM indexes are supported. Since there is not a full CRAM parser available in go yet, `indexcov` uses only
the .crai files and requires a `.fasta.fai` to be sent via `--fai` so that it knows the chromosome names around
lengths. The sample names are inferred from the file names. crai resolution is often much less than 100KB (compared to)
16KB for the bam index, but it is sufficient to find large-scale differences in coverage.

Example usage with cram looks like:

```
goleft indexcov -d output/ --fai h human_g1k_v37.fasta.fai /path/to/*.crai
```

**note** that the .fai (not the fasta) is required and that the files are .crai (not cram).

How It Works
============

The bam index stores a linear index for each chromosome indicating the file (and virtual) offset for every 16,384 bases in
that chromosome. Since we know the total number of 16,384 base intervals in the index and the size of the bam file (from the
last file offset stored in the index), we know the average size (in bytes) of taken by each 16,384 base chunk. So, we iterate
over each (16KB) element in the linear index, subtract the previous file offset, and scale by the expected (average) size. This
gives the scaled value for each 16,384-base chunk. There are many ways that this value can be off, but, in practice, it works
well as a rough estimate.

Because of this `indexcov` is of less-use on exome or targetted capture, but those will
be very fast to run with `goleft depth` anyway.

<a name="Files"></a> Files
==========================

In addition to the  interactive HTML files, `indexcov` outputs a number of text files:

+ `$prefix-indexcov.ped`: a .ped/.fam file with the inferred sex in the appropriate column if the sex chromosomes were found.
                          the CNX and CNY columns indicating the floating-point estimate of copy-number for those chromosomes.
                          `bins.out`: how many bins had a coverage value outside of (0.85, 1.15). high values can indicate high-bias samples.
                          `bins.lo`: number of bins with value < 0.15. high values indicate missing data.
                          `bins.hi`: number of bins with value > 1.15. 
                          `bins.in`: number of bins with value inside of (0.85, 1.15)
                          `p.out`: `bins.out/bins.in`
                          `PC1...PC5`: PCA projections calculated with depth of autosomes.

+ `$prefix-indexcov.roc`: tab-delimited columns of chrom, scaled coverage cutoff, and $n_samples columns where each indicates the
                          proportion of 16KB blocks at or above that scaled coverage value.
+ `$prefix-indexcov.bed.gz`: a bed file with columns of chrom, start, end, and a column per sample where the values indicate there
                             scaled coverage for that sample in that 16KB chunk.
