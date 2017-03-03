indexsplit
----------

`indexsplit` quickly generates evenly sized (by amount of data) regions across
a cohort. It does this by reading the bam (or cram) index and using the file offsets as proxies
for the amount of data. It sums the values in these bins across all samples. This gives a good
estimate for actual reads in the region but without having to parse the bam file.

A common use of this will be to generate regions to be used to parallelize variant calling fairly
by splitting in to `N` regions with approximately equal amounts of data **across the cohort**.


Usage
-----

```
goleft indexsplit -N 5000 --fai hs37d5.fa.fai /path/to/*.bam > regions.bed
```

TODO
----

+ evaluate sequence complexity of each region (using fasta) and split low-complexity regions further.
