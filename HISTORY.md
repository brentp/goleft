v0.1.10 (dev)
=============

+ report correct version # (v0.1.9 reported 0.1.8 from CLI)
+ `indexcov`: bigger hit radius in plots == easier to mouse-over points
+ `indexcov`: use vOffset for better resolution of (estimated) depth.
+ `indexcov`: ~2X faster due to optimizations in biogo/hts and buffered IO.

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
