indexcov
========

quickly estimate coverage from the bam index. 
A bam index has 16KB resolution so that's what this gives, but it gives what appears to be a high-quality estimate
in seconds per genome.

The output is scaled to around 1. So a long stretch with values of 2 would be a duplication.
This is useful as a quick QC to get coverage values across the genome.

Usage
=====

```
goleft indexcov -c $chrom *.bam > depth.bed
```

This will create a bed file where each additional column is the normalized, estimated depth for each
sample.

With that, it's simple to plot the depth across the genome:

![Example](https://cloud.githubusercontent.com/assets/1739/21273832/a42c3a6c-c382-11e6-9bd1-3870a8333c04.png "example depth")

Where each color is a sample and here we can see that the purple sample has a large deletion and the yellow/tan sample has a
higher variance.


