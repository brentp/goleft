#!/bin/bash

goleft depth --windowsize 10 --q 1 --mincov 4 --reference hg19.fa --processes 1 --stats --bed Test1-coverage.depth-tocalculate-windows.bed --prefix out Test1-sort.bam
