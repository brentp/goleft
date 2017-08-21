ref=Homo_sapiens.GRCh38.dna.toplevel.fa
bam=NA12878_S1.bam
set -eu
chrom=chr1
if [[ ! -e st.depth.$chrom.depth.bed ]]; then
time ~/bin/goleft depth --chrom $chrom -p 20 -o -w 16384 --reference $ref $bam --prefix st.depth
#time samtools depth $bam > /dev/null
fi
time ~/bin/goleft indexcov -d ix.depth --sex "chrX,chrY" $bam
zgrep -w $chrom ix.depth/ix.depth-indexcov.bed.gz > ix.depth/$chrom.bed
python cmp.py ix.depth/$chrom.bed st.depth.$chrom.depth.bed


bedtools intersect -sorted -a st.depth.$chrom.depth.bed -b <(zgrep -w "^chr1" LCR-hs38.bed.gz) -c > noLCR.st.bed
bedtools intersect -sorted -a ix.depth/$chrom.bed -b <(zgrep -w "^chr1" LCR-hs38.bed.gz) -c > noLCR.ix.bed
python cmp.py noLCR.ix.bed noLCR.st.bed LCR
