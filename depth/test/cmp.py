import subprocess as sp
import sys

goleft_bed = sys.argv[1]
bam = sys.argv[2]

for toks in (l.rstrip().split("\t") for l in open(goleft_bed)):
    out = sp.check_output("samtools depth -Q 1 -r '%s:%d-%s' %s | awk '{s+=$3}END{if(NR==0){print 0}else{print s/NR}}'" % (toks[0], int(toks[1]) + 1, toks[2], bam), shell=True).strip()
    expected = float(toks[3])

    if expected - float(out.strip()) > 0.5:
        print("ERROR")
        print(expected - float(out.strip()), expected)
        print("samtools depth -Q 1 -r '%s:%d-%s' %s | awk '{s+=$3}END{print s/NR}'" % (toks[0], int(toks[1]) - 1, toks[2], bam))
        sys.exit(1)

