Bins
====

A "*bin*" in `indexcov` terms is a 16,384 base window for which coverage was estimated.
For a *good* sample, we expect most of the bins to be around 1.

This plot shows:
+ x-axis: number of bins with scaled coverage < 0.15
+ y-axis: number of bins with scaled coverage outside of 0.85 - 1.15

These relatively simple metrics are quite effective in finding problematic samples.

samples that fall:
+ to the left of the plot have many regions with low or missing coverage.
+ high on the plot have many regions outside of the expected coverage interval.

The former of those may be due to truncated bam files or other missing data.
The latter is a good indicator of dosage bias (basically uneven coverage across the genome).

Samples that have large values on the y-axis will be problematic for CNV and SV calling.

For example:
![bin Help](https://cloud.githubusercontent.com/assets/1739/22121227/6ccb2ffa-de40-11e6-8916-eb7f3a584c35.png "bin help")

In this plot, we see 2 samples that have very extreme values for both the x and y axis. These are samples
that were missing coverage for entire sections of the genome due to processing glitches.

There are also several samples with values above 25000 on the y-axis that warrant further investigation.

We can not use hard cutoffs like 25000 because the number of regions covered will vary according to many factors
such as the reference genome used, the mapper, the PCR artefacts, etc. It is best to evaluate based on the
cohort.
