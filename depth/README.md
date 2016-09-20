depth
=====

depth parallelizes calls to [samtools](https://samtools.github.io) in user-defined windows.
It outputs a bed file of callable regions (determined by mincov) and of depth (only windows
with >= `maxmeandepth` are reported.

```
usage: goleft depth [--windowsize WINDOWSIZE] [--maxmeandepth MAXMEANDEPTH] [--q Q] [--chrom CHROM] [--mincov MINCOV] [--stats] --reference REFERENCE [--processes PROCESSES] [--bed BED] [--prefix PREFIX] BAM

positional arguments:
  bam                    bam for which to calculate depth

options:
  --windowsize WINDOWSIZE, -w WINDOWSIZE
                         window size in which to calculate high-depth regions [default: 250]
  --maxmeandepth MAXMEANDEPTH, -m MAXMEANDEPTH
                         windows with depth > than this are high-depth. The default reports the depth of all regions.
  --q Q, -Q Q            mapping quality cutoff [default: 1]
  --chrom CHROM, -c CHROM
                         optional chromosome to limit analysis
  --mincov MINCOV        minimum depth considered callable [default: 4]
  --stats, -s            report sequence stats [GC CpG masked] for each window
  --reference REFERENCE, -r REFERENCE
                         path to reference fasta
  --processes PROCESSES, -p PROCESSES
                         number of processors to parallelize.
  --bed BED, -b BED      file of positions or regions.
  --prefix PREFIX
  --help, -h             display this help and exit
