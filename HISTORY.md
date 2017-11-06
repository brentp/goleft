v0.1.18 (dev)
=======
+ `indexsplit`: fix off-by-one that resulted in double-counting some regions.
+ `indexcov`: cram edge-cases.
+ `indexcov`: better normalization to 1 for all cases. Fixes bug for bams with
              many (e.g. > 10K) chromosomes of which many have very low or normalization
              coverage. (#36)

v0.1.17
=======
+ `indexcov`: dont error when no sex chromosomes are found (#27).
+ `indexcov`: dont error when some chromosomes have a single region and others have 0.
+ `indexcov`: better checking on short sex chroms and other CRAI fixes.

v0.1.16
=======
+ `indexcov`: report and plot number of mapped and unmapped reads as reported by the index.
+ `covmed`: rename to `covstats`
+ `covstats`: report samplename(s) from read groups as well as bam path.
+ `covstats`: skip first 100K reads to give better estimates of depth.
+ `covstats`: report percent of bad (QC-Fail|Duplicate) and of umapped reads.
+ `indexcov`: automatically exclude chromosomes that match pattern: `^chrEBV$|^NC|_random$|Un_|^HLA\-|_alt$|hap\d$`.
              this can be adjust from the command-line.

v0.1.15
=======
+ `indexsplit`: more closely match request number of regions.
+ `indexsplit`: allow specifying "problematic" regions which are split more finely.
+ `samplename`: **new tool** report the sample name of a bam from the read-group in its header.

v0.1.14
=======
+ `indexcov`: speed and memory improvements.
+ `indexcov`: support crai.
+ `indexcov`: prettier plot and labels for pngs.
+ `indexcov`: bin plot shows proportions rather than total numbers.
+ `indexcov`: better messages and docs for `--sex` argument.
+ `indexsplit`: new tool to generate regions for parallelization (from bam/cram indexes)

v0.1.13
=======
+ `covmed`: better coverage estimate by scaling by dup|qcfail|secondary.
+ `covmed`: report 5th and 95th percentile of insert size.
+ `indexcov`: make index.html even when no sex chromosomes are present.
+ `depth:` quote reference and chromosome name.

v0.1.12
=======
+ `indexcov`: --sex is now specified as a comma-delimited string
+ `indexcov`: handle a single sample (see #16)
+ `indexcov`: functional tests.
+ `indexcov`: add "slope" output which indicates the slope of the coverage plot between ~0.85 and ~1.15.
              this matches bins.out fairly well, but it is another metric to look at.

v0.1.11
=======
+ `indexcov`: make plots a saner size.
+ `indexcov`: output an index.html page that summarizes findings and links to other pages.
              the index contains static pngs for coverage and depth that link to the interactive
              HTML pages for each chromosome.
+ `indexcov`: remove --prefix argument in favor of -d/--directory and write the overview page index
              index.html


v0.1.10
=======

+ report correct version # (v0.1.9 reported 0.1.8 from CLI)
+ `indexcov`: bigger hit radius in plots == easier to mouse-over points
+ `indexcov`: use vOffset for better resolution of (estimated) depth.
+ `indexcov`: ~2X faster due to optimizations in biogo/hts and buffered IO.
+ `indexcov`: PCA and dosage bias plots (values also reported in .ped file).

v0.1.9
======
+ `indexcov`: fix for sparse bams and other edge-cases.
+ `indexcov`: performance improvements.
+ `indexcov`: output html with plot of coverage.
+ `indexcov`: plot and report sex inferred from coverage on allosomes.

v0.1.8
======
+ `indexcov`: output header with sample names.
+ `indexcov`: limit reported CN.

v0.1.7
======

+ fix bug in goleft depth where some lines were output multiple times.
+ add covmed subprogram to calculate median coverage, mean and sd of insert size quickly.
+ add depthwed to make a matrix file from multiple .depth files from `goleft depth`.
+ goleft depth will always output windows of the exact size requested, previously
  if 10MB was not a multiple of windowsize, we'd get smaller windows at the edge.

+ new tool: **indexcov**: very fast calculation of relative coverage using the bam index.
