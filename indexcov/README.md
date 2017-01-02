indexcov
========

Quickly estimate coverage from the bam index. 
A bam index has 16KB resolution so that's what this gives, but it provides what appears to be a high-quality 
coverage estimate in seconds per genome.

The output is scaled to around 1. So a long stretch with values of 1.5 would be a heterozygous duplication.
This is useful as a quick QC to get coverage values across the genome.

In our tests, we can **estimate depth across a 60X genomes for 30 samples in 1 minute**.

Interactive HTML plots of depth are output for each chromosome.

Usage
=====

```
goleft indexcov -c $chrom *.bam > depth.bed
```

This will create a bed file where each additional column is the normalized, estimated depth for each
sample. If no chromosome is given, it will do the whole genome. Since most of the time is scaling There
coverages based on the entire index, there's very little speed benefit to choosing only a single chromosome.

In addition, an interactive HTML file is output per chromosome that displays information like this:

![Example](https://cloud.githubusercontent.com/assets/1739/21273832/a42c3a6c-c382-11e6-9bd1-3870a8333c04.png "example depth")

Where each color is a sample and here we can see that the purple sample has a large deletion and the yellow/tan sample has a
higher variance.

If we view the html file for **X chromosome** we can see a nice separation of samples by sex along with the PAR at the left:

![X Example](https://cloud.githubusercontent.com/assets/1739/21597648/074f06ca-d10b-11e6-8732-e9a2e8d1ecb5.png "sex example")

That plot is taken directly from the HTML output by `indexcov`.

Finally, `indexcov` will output a coverage (ROC) plot that shows how much of the genome is coverage at at given (scaled) depth.
This is output to a $prefix-depth-roc.html file and looks like:
![ROC Example](https://cloud.githubusercontent.com/assets/1739/21597688/ef76ae8a-d10b-11e6-88a0-735c73203fc7.png "ROC example")

Here we can see that one sample has much lower coverage than the rest, and we can hover and determine the exact sample.


How It Works
============

The bam index stores a linear index for each chromosome indicating the file (and virtual) offset for every 16,384 bases in
that chromosome. Since we know the total number of 16,384 base intervals in the index and the size of the bam file (from the
last file offset stored in the index), we know the average size (in bytes) of taken by each 16,384 base chunk. So, we iterate
over each (16KB) element in the linear index, subtract the previous file offset, and scale by the expected (average) size. This
gives the scaled value for each 16,384-base chunk. There are many ways that this value can be off, but, in practice, it works
well as a rough estimate.

Todo
====

parallelize the depth calc by sample.
