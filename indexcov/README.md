indexcov
========

quickly estimate coverage from the bam index. 
A bam index has 16KB resolution so that's what this gives, but it gives what appears to be a high-quality estimate
in seconds per genome.

The output is scaled to around 1. So a long stretch with values of 2 would be a duplication.

Usage
=====

```
goleft indexcov $chrom $bam > depth.bed
```
