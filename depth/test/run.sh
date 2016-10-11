#!/bin/bash
set -eu -o pipefail

# Standard run
goleft depth --windowsize 10 --q 1 --mincov 4 --reference hg19.fa --processes 1 --stats --bed windows.bed --prefix out t.bam


# Empty BAM file, simulating empty chromosomes
goleft depth --windowsize 10 --q 1 --mincov 4 --reference hg19.fa --processes 1 --stats --prefix out-empty t-empty.bam
