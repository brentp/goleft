#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest
set -uo pipefail

go build -o ./goleft_test ../cmd/goleft/goleft.go

run check_sample_name ./goleft_test samplename ../indexcov/samples/sample_paper_0021.bam
assert_exit_code 0
assert_in_stdout sample_paper_0021

run check_sample_name_e ./goleft_test samplename -e ../indexcov/samples/sample_paper_0021.bam
assert_exit_code 0
assert_in_stdout sample_paper_0021

