import subprocess as sp
import sys

goleft_bed = sys.argv[1]
bam = sys.argv[2]

for toks in (l.rstrip().split("\t") for l in open(goleft_bed)):
    cmd = "samtools depth -a -Q 1 -r '%s:%d-%s' %s | awk '{s+=$3}END{if(NR==0){print 0}else{print s/%d}}'" % (toks[0], int(toks[1]) + 1, toks[2], bam, int(toks[2]) - int(toks[1]))
    out = sp.check_output(cmd, shell=True).strip()
    expected = float(toks[3])

    if abs(expected - float(out.strip())) > 0.5:
        print("ERROR")
        print(float(out.strip()), expected)
        print(cmd)
        sys.exit(1)

