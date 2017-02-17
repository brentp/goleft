import sys

freal = sys.argv[1]
fcram = sys.argv[2]

def read(path):
    xs, ys = [], []
    for toks in (l.rstrip().split("\t") for l in open(path)):
        try:

            ys.append(float(toks[3]))
            xs.append(int(toks[1]))
        except IndexError:
            break
    return xs, ys 

realx, realy = read(freal)
cramx, cramy = read(fcram)

n = min(len(realx), len(cramx))

import numpy as np

realx, realy = realx[:n], np.array(realy[:n])
cramx, cramy = cramx[:n], np.array(cramy[:n])


cramy = cramy.astype(float) / np.median(np.array(cramy))
realy = realy.astype(float) / np.median(realy)

assert all(c == r for c, r in zip(cramx, realx))

from matplotlib import pyplot as plt

import seaborn as sns
sns.set_palette('Set1')
sns.set_style('whitegrid')

plt.plot(realx, realy, marker=None, lw=1, ls='-', label='samtools-depth')
plt.plot(realx, cramy, marker=None, lw=1, ls='-', label='cram-index')
plt.ylim(0, 5)
plt.xlabel("genomic position")
plt.ylabel("scaled depth")
plt.legend()
plt.show()
