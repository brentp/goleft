import sys
from itertools import groupby
import numpy as np


def gen():
    for l in sys.stdin:
        chrom, p, d = l.rstrip().split()
        if chrom != '1': break
        p, d = int(p), float(d)
        yield p, d

last = 0
interval = 16384
for i, g in groupby(gen(), lambda (p, d) : p / interval):
    if i * interval > 150028288:
        break

    g = list(g)
    n = len(g)
    vals = [x[1] for x in g]
    while last < i * interval:
        print "1\t%d\t%d\t0" % (last, last+16384)
        last += 16384

    for k in range(interval/16384):
        print "1\t%d\t%d\t%.2f" % ((i+k) * 16384, (i+k+1)*16384, sum(vals) / 16384)
    last = (i+1) * interval
