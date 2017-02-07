import sys
import gzip
import itertools as it
import numpy as np
import scipy.stats as ss
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

fha = (gzip.open if sys.argv[1].endswith(".gz") else open)(sys.argv[1])
fhb = (gzip.open if sys.argv[2].endswith(".gz") else open)(sys.argv[2])

def gen(fh):
    for line in fh:
        toks = line.rstrip().split("\t")
        toks[1], toks[2] = int(toks[1]), int(toks[2])
        toks[3] = float(toks[3])
        yield toks

xs, ys = [], []

for i, (a, b) in enumerate(it.izip(gen(fha), gen(fhb))):

    if a[1] != b[1]:
        raise Exception("expected same positions for both files")

    xs.append(a[3])
    ys.append(b[3])
    if xs[-1] < 0.05 and ys[-1] > 20:
        print(a, b)

fig, axes = plt.subplots(1)
axes = (axes,)

ys = np.array(ys)
ys /= np.median(ys)
xs = np.array(xs)

diff = xs - ys
out = sum(abs(d) > 0.5 for d in diff)
print "out:", out, "total:", len(diff),  ("%.2f" % (100.0*out/len(diff)))

print sum(abs(d) < 0.25 for d in diff) / float(len(diff))
diff = diff[np.abs(diff) < 0.5]

axes[0].hist(diff, 20)
axes[0].set_xlim(-0.5, 0.5)
axes[0].set_xlabel("Difference in depth estimate (indexcov - samtools)")
axes[0].set_ylabel("Count")

d = "/uufs/chpc.utah.edu/common/home/u6000771/public_html/"
plt.savefig(d + "t.png")
plt.show()
