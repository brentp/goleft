v0.1.7
======

+ fix bug in goleft depth where some lines were output multiple times.
+ add covmed subprogram to calculate median coverage, mean and sd of insert size quickly.
+ add depthwed to make a matrix file from multiple .depth files from `goleft depth`.
+ goleft depth will always output windows of the exact size requested, previously
  if 10MB was not a multiple of windowsize, we'd get smaller windows at the edge.

+ new tool: **indexcov**: very fast calculation of relative coverage using the bam index.
