indexsplit
----------

`indexsplit` quickly generates evenly sized (by amount of data) regions across
a cohort. It does this by reading the bam (or cram) index and using the file offsets as proxies
for the amount of data. It sums the values in these bins across all samples. This gives a good
estimate for actual reads in the region but without having to parse the bam file.

The result is a bed file with an additional column indicating the (scaled) size of data in these
region across all samples. The numbers in that column will be fairly even except at the ends offsets
chromosomes or for small chromosomes.

A common use of this will be to **generate regions to be used to parallelize variant-calling** fairly
by splitting in to `N` regions with approximately equal amounts of data **across the cohort**.

On a modest laptop with an SSD, `indexsplit` can **generate even-coverage regions in ~4 seconds for 45 bams**.
The time is independent of the number of regions.

When a single 16KB chunk will has more data than the determined chunk size, `indexpslit` will output
sub-regions of that chunk even though it doesn't know the exact placement of the data within it.

Usage
-----

```
goleft indexsplit -N 5000 /path/to/*.bam > regions.bed
```

If you want to do this from CRAM files send the .crai files and a fasta index via
the `--fai`:

```
goleft indexsplit -N 8000 --fai reference.fa.fai /path/to/*.crai > regions.bed
```

The user is responsible for unsuring that the crai chromosome order matches the .fai order 
(this will be the case if the fasta was the same as used in alignment).
