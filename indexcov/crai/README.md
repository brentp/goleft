
read a crai and output sizes.

the python scripts in this directory simplify comparing the output of samtools depth

with the output of crai.

e.g.

```
samtools depth xx/phase3/data/NA21144/alignment/NA21144.mapped.ILLUMINA.bwa.GIH.low_coverage.20130415.bam.cram > t.out
cat t.out | python agg.py > real.bed
python plot-vs-with-st-depth.py real.bed o.bed
```
where o.bed is from the crai.Sizes() (and the output of running go test in this directory).
