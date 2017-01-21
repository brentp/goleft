#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest
set -uo pipefail

go build -o ./goleft_test ../cmd/goleft/goleft.go

check_nobams() {
    echo "XXX"
}


export -f check_nobams

run check_no_bams ./goleft_test indexcov -d /tmp/tt
assert_exit_code 255
assert_in_stderr "expected at least 1 bam"
