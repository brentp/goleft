#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest
set -e

go build -o goleft ../main.go

check_with_fai_bt() {
    fai=$1
    bed=$2
    v=$(awk '{ print $1"\t"0"\t"$2 }' $fai | bedtools subtract -b - -a $bed && awk '{ print $1"\t"0"\t"$2 }' $fai | bedtools subtract -a - -b $bed)
    echo -e "$v"
}

check_with_bed_bt() {
    v=$(bedtools subtract -a $1 -b $2 && bedtools subtract -a $2 -b $1)
    echo -e "$v"
}
export -f check_with_fai_bt


run check_wgs ./goleft depth -Q 1 --ordered --windowsize 10 --stats --prefix x --reference test/hg19.fa test/t.bam

assert_exit_code 0
assert_equal "$(check_with_fai_bt test/hg19.fa.fai x.depth.bed)" ""
assert_equal "$(check_with_fai_bt test/hg19.fa.fai x.callable.bed)" ""


run check_wgs_big_window ./goleft depth -Q 1 --ordered --windowsize 1000000000 --stats --prefix x --reference test/hg19.fa test/t.bam
assert_exit_code 0
assert_equal "$(check_with_fai_bt test/hg19.fa.fai x.depth.bed)" ""
assert_equal "$(check_with_fai_bt test/hg19.fa.fai x.callable.bed)" ""

for w in 55 60 71 13 2001; do
    run check_wgs_big_window$w ./goleft depth -Q 1 --ordered --windowsize $w --stats --prefix x --reference test/hg19.fa test/t.bam
    assert_exit_code 0
    assert_equal "$(check_with_fai_bt test/hg19.fa.fai x.depth.bed)" ""
    assert_equal "$(check_with_fai_bt test/hg19.fa.fai x.callable.bed)" ""
done

run check_bed ./goleft depth --bed test/windows.bed -Q 1 --ordered --windowsize 10 --stats --prefix x --reference test/hg19.fa test/t.bam
assert_exit_code 0
assert_equal "$(check_with_bed_bt x.depth.bed test/windows.bed)" ""
assert_equal "$(check_with_bed_bt x.callable.bed test/windows.bed)" ""

run check_bed_big_window ./goleft depth --bed test/windows.bed -Q 1 --ordered --windowsize 1000000 --stats --prefix x --reference test/hg19.fa test/t.bam
assert_exit_code 0
assert_equal "$(check_with_bed_bt x.depth.bed test/windows.bed)" ""
assert_equal "$(check_with_bed_bt x.callable.bed test/windows.bed)" ""

for w in 50 55 60 71 13 2002; do
    run check_bed_window$w ./goleft depth --bed test/windows.bed -Q 1 --ordered --windowsize $w --stats --prefix x --reference test/hg19.fa test/t.bam
    assert_exit_code 0
    assert_equal "$(check_with_bed_bt x.depth.bed test/windows.bed)" ""
    assert_equal "$(check_with_bed_bt x.callable.bed test/windows.bed)" ""
done
