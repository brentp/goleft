#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest
set -uo pipefail

go build -o ./goleft_test ../cmd/goleft/goleft.go

check_nobams() {
    echo "XXX"
}

get() {
    set -e
    bam=$1
    if [[ ! -e $bam ]]; then
        curl -sS https://s3.amazonaws.com/b4-test-data/$bam > $bam
    fi
    if [[ ! -e $bam.bai ]]; then
        curl -sS https://s3.amazonaws.com/b4-test-data/$bam.bai > $bam.bai
    fi
    set +e
}


export -f check_nobams
run check_no_bams ./goleft_test indexcov -d /tmp/tt
assert_exit_code 255
assert_in_stderr "expected at least 1 bam"


get sample_name_0001.bam
rm -rf /tmp/tt
run check_single_sample ./goleft_test indexcov -d /tmp/tt sample_name_0001.bam
assert_exit_code 0
assert_in_stderr "not plotting"

run check_sex ./goleft_test indexcov -d /tmp/tt --sex "X,Y" sample_name_0001.bam
assert_exit_code 0

get sample_biogo_e_0001.bam
run check_issue17 ./goleft_test indexcov -d /tmp/tt sample_biogo_e_0001.bam
assert_exit_code 1
assert_in_stderr "no usable chromsomes in bam"

run check_sex_warning ./goleft_test indexcov --sex chrX,chrY -d /tmp/tt sample_name_0001.bam
assert_exit_code 1
assert_in_stderr "found chromosome \"X\", wanted \"chrX\""
assert_in_stderr "found chromosome \"Y\", wanted \"chrY\""

mkdir -p samples/ && cd samples/
if [[ ! -e sample_paper.tar.xz ]]; then
    curl -Ss https://s3.amazonaws.com/b4-test-data/sample_paper.tar.xz > sample_paper.tar.xz
fi
if [[ ! -e sample_paper_0030.bam.bai ]]; then
    tar xJf sample_paper.tar.xz
fi
cd ../

run check_cohort ./goleft_test indexcov -d /tmp/tt samples/*.bam
assert_exit_code 0


rm -r /tmp/tt
run check_no_sex ./goleft_test indexcov --sex "" -d /tmp/tt samples/*.bam
assert_in_stderr "index.html for overview"

rm samples/*.bam
