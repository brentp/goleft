from matplotlib import pyplot as plt
plt.rcParams['pdf.fonttype'] = 42

import pandas as pd
import seaborn as sns
sns.set_style('whitegrid')
sns.set_palette('Set1')

pcr_plus = set(x.strip() for x in open('pcr-plus.samples'))

df = pd.read_table('simons-indexcov.ped', index_col='sample_id')
size = 9

pp = df.loc[pcr_plus]

tot = pp['bins.in'] + pp['bins.out']

plt.scatter(pp['bins.lo'] / tot, pp['bins.out'] / tot, alpha=0.55, label='PCR +', s=size, zorder=2)

pf = df.loc[[x for x in df.index if not x in pcr_plus]]
tot = pf['bins.in'] + pf['bins.out']
plt.scatter(pf['bins.lo'] / tot, pf['bins.out'] / tot, alpha=0.55, label='PCR free', s=size)

plt.xlabel('Proportion of bins with scaled coverage < 0.15')
plt.ylabel('Proportion of bins with scaled coverage outside of (0.85, 1.15)')
plt.xlim(xmin=0.067, xmax=0.0695)

plt.legend()
plt.savefig('figure4.eps', dpi=1000)
