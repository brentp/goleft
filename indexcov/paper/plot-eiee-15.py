from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
sns.set_style('white')

df = pd.read_table('eiee.15.bed.gz', compression='gzip')
print(df.head())
cols = list(df.columns)
cols[0] = "chrom"
cols = [c for c in cols[3:] if c != '15-0022964']


fig, ax = plt.subplots(1)

for c in cols:
    ax.plot(df['start'], df[c], color='#dddddd', lw=0.4)

ax.plot(df['start'], df['15-0022964'], color='#487535', lw=0.4)
ax.set_ylabel("Scaled Coverage")
ax.set_xlabel("Position On Chromosome 15")
ax.set_xlim(xmin=0, xmax=df.start.max())
print(df.start.max())
ax.set_ylim(0, 6)
plt.draw()
ticks = ax.get_xticks()

labels = ["%dM" % (t / 1000000) for t in ticks if t < df.start.max()]
ax.set_xticks(ticks, labels)
ax.set_xticklabels(labels)
ax.set_xlim(xmin=0, xmax=df.start.max())
sns.despine()

plt.show()

plt.close()


fig, ax = plt.subplots(1)

df = pd.read_table('eiee.15.roc')
for c in cols:
    ax.plot(df['cov'], df[c], color='#dddddd', lw=2)

ax.plot(df['cov'], df['15-0022964'], color='#487535', lw=2)
ax.set_xlabel("Scaled Coverage")
ax.set_ylabel("Proportion of Regions Covered")
sns.despine()
plt.show()
