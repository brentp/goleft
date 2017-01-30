#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest
set -eo pipefail

go build -o goleft ../cmd/goleft/goleft.go

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

check_uniq() {
    a=$(cat $1 | wc -l)
    b=$(uniq $1 | wc -l)
    if [[ "$a" != "$b" ]]; then
        echo "DUPLICATE LINES in $1: $a $b"
        return
    fi
    if [[ "$#" == "2" ]]; then
        echo "OK"
        return
    fi
    c=$(cut -f 1-3 $1 | sort -u | wc -l)
    if [[ "$a" != "$c" ]]; then
        echo "DUPLICATE regions in $1: $a $c"
        return
    fi
    echo "OK"
}

export -f check_with_fai_bt
export -f check_with_bed_bt


run check_wgs ./goleft depth -Q 1 --ordered --windowsize 100 --stats --prefix x --reference test/hg19.fa test/t.bam

assert_exit_code 0
assert_equal "$(check_with_fai_bt test/hg19.fa.fai x.depth.bed)" ""
assert_equal "$(check_with_fai_bt test/hg19.fa.fai x.callable.bed)" ""
assert_equal "$(check_uniq x.depth.bed)" "OK"
assert_equal "$(check_uniq x.callable.bed)" "OK"

run compare_to_samtools_100 python test/cmp.py x.depth.bed test/t.bam
assert_exit_code 0

run check_wgs_big_window ./goleft depth -Q 1 --ordered --windowsize 1000000000 --stats --prefix x --reference test/hg19.fa test/t.bam
assert_exit_code 0
assert_equal "$(check_with_fai_bt test/hg19.fa.fai x.depth.bed)" ""
assert_equal "$(check_with_fai_bt test/hg19.fa.fai x.callable.bed)" ""
assert_equal "$(check_uniq x.depth.bed)" "OK"
assert_equal "$(check_uniq x.callable.bed)" "OK"

for w in 55 60 71 13 2001; do
    run check_wgs_big_window$w ./goleft depth -Q 1 --ordered --windowsize $w --stats --prefix x --reference test/hg19.fa test/t.bam
    assert_exit_code 0
    assert_equal "$(check_with_fai_bt test/hg19.fa.fai x.depth.bed)" ""
    assert_equal "$(check_with_fai_bt test/hg19.fa.fai x.callable.bed)" ""
    assert_equal "$(check_uniq x.depth.bed)" "OK"
    assert_equal "$(check_uniq x.callable.bed)" "OK"
done


run check_bed ./goleft depth --bed test/windows.bed -Q 1 --ordered --windowsize 10 --stats --prefix x --reference test/hg19.fa test/t.bam
assert_exit_code 0
assert_equal "$(check_with_bed_bt x.depth.bed test/windows.bed)" ""
assert_equal "$(check_with_bed_bt x.callable.bed test/windows.bed)" ""
assert_equal "$(check_uniq x.depth.bed bed)" "OK"
assert_equal "$(check_uniq x.callable.bed bed)" "OK"


run check_bed_big_window ./goleft depth --bed test/windows.bed -Q 1 --ordered --windowsize 1000000 --stats --prefix x --reference test/hg19.fa test/t.bam
assert_exit_code 0
assert_equal "$(check_with_bed_bt x.depth.bed test/windows.bed)" ""
assert_equal "$(check_with_bed_bt x.callable.bed test/windows.bed)" ""
assert_equal "$(check_uniq x.depth.bed bed)" "OK"
assert_equal "$(check_uniq x.callable.bed bed)" "OK"
run compare_to_samtools_big_window python test/cmp.py x.depth.bed test/t.bam

for w in 50 55 60 71 13 2002; do
    run check_bed_window$w ./goleft depth --bed test/windows.bed -Q 1 --ordered --windowsize $w --stats --prefix x --reference test/hg19.fa test/t.bam
    assert_exit_code 0
    assert_equal "$(check_with_bed_bt x.depth.bed test/windows.bed)" ""
    assert_equal "$(check_with_bed_bt x.callable.bed test/windows.bed)" ""
    assert_equal "$(check_uniq x.depth.bed bed)" "OK"
    assert_equal "$(check_uniq x.callable.bed bed)" "OK"
    run compare_to_samtools_directly$W python test/cmp.py x.depth.bed test/t.bam
done

assert_exit_code 0


run check_empty ./goleft depth --windowsize 10 --q 1 --mincov 4 --reference test/hg19.fa --processes 1 --stats --prefix x test/t-empty.bam
assert_exit_code 0
assert_equal "$(check_with_fai_bt test/hg19.fa.fai x.depth.bed)" ""
assert_equal "$(check_with_fai_bt test/hg19.fa.fai x.callable.bed)" ""
assert_equal "$(check_uniq x.depth.bed)" "OK"
assert_equal "$(check_uniq x.callable.bed)" "OK"


run check_empty_window ./goleft depth --bed test/windows.bed --windowsize 10 --q 1 --mincov 4 --reference test/hg19.fa --processes 1 --stats --prefix x test/t-empty.bam
assert_exit_code 0
assert_equal "$(check_with_bed_bt x.depth.bed test/windows.bed)" ""
assert_equal "$(check_with_bed_bt x.callable.bed test/windows.bed)" ""
assert_equal "$(check_uniq x.depth.bed bed)" "OK"
assert_equal "$(check_uniq x.callable.bed bed)" "OK"


run check_hla ./goleft depth -r test/fake.fa --prefix /tmp/xx test/hla.bam
assert_exit_code 0

echo -e "\nFINISHED OK"

