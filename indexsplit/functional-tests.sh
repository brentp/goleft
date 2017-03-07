#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest
set -uo pipefail

go build -o ./goleft_test ../cmd/goleft/goleft.go

mkdir -p ../indexcov/samples/ && cd ../indexcov/samples/
if [[ ! -e sample_paper.tar.xz ]]; then
    curl -Ss https://s3.amazonaws.com/b4-test-data/sample_paper.tar.xz > sample_paper.tar.xz
fi
if [[ ! -e sample_paper_0030.bam.bai ]]; then
    tar xJf sample_paper.tar.xz
fi
if [[ ! -e sample_paper_0030.bam ]]; then
    tar xJf sample_paper.tar.xz
fi
cd -

if [[ ! -e human_g1k_v37.fasta.fai ]]; then
    wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
fi

splitter() {
    ./goleft_test indexsplit -n 1000 --fai human_g1k_v37.fasta.fai ../indexcov/samples/*.bam > _test.cov
}
samtools view -H ../indexcov/samples/sample_paper_0001.bam | grep "^@SQ" | perl -pe 's/.+SN://' | awk 'BEGIN{FS="\tLN:"; OFS="\t"} { print $1,0,$2 }' > genome.bed

run check_it splitter
assert_exit_code 0
cat $STDERR_FILE
assert_equal 0 $(bedtools subtract -a genome.bed -b _test.cov | wc -l)
assert_equal 0 $(bedtools subtract -b genome.bed -a _test.cov | wc -l)
assert_equal 0 $(bedtools intersect -a _test.cov -b _test.cov -c | awk '$NF != 1' | wc -l)
rm -f genome.bed _test.cov
