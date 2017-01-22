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

