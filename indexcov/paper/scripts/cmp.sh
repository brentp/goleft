ref=Homo_sapiens.GRCh38.dna.toplevel.fa
bam=NA12878_S1.bam
set -eu
chrom=chr1
if [[ ! -e st.depth.$chrom.depth.bed ]]; then
time ~/bin/goleft depth -p 20 --chrom $chrom -o -w 16384 --reference $ref $bam --prefix st.depth
fi
time ~/bin/goleft indexcov -d ix.depth --sex "chrX,chrY" $bam
zgrep -w $chrom ix.depth/ix.depth-indexcov.bed.gz > ix.depth/$chrom.bed
python cmp.py ix.depth/$chrom.bed st.depth.$chrom.depth.bed
