from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
sns.set_style('whitegrid')

df = pd.read_table('simons-indexcov.with-pca-free.ped')

fig, ax = plt.subplots(1)
"""
family_id	sample_id	paternal_id	maternal_id	sex	phenotype	CNX	CNY	bins.out	bins.lo	bins.hi	bins.in	slope	p.out	PC1	PC2	PC3	PC4	PC5	pcr
unknown	SS0012978	-9	-9	2	-9	1.99	0.00	16893	11865	3515	159270	118.575	0.11	-206.56	404.37	1830.56	-864.83	2393.05	0
unknown	SS0012979	-9	-9	1	-9	1.02	1.03	16803	11830	3344	159360	118.681	0.11	1673.99	-705.83	2592.70	-2354.44	-624.61	0
unknown	SS0013012	-9	-9	2	-9	1.99	0.00	17248	11848	3716	158915	118.316	0.11	-965.45	-867.86	1960.27	-1007.85	1912.20	0
unknown	SS0013018	-9	-9	1	-9	1.03	1.03	16884	11840	3393	159279	118.592	0.11	1132.93	-1922.86	1063.57	-2462.44	590.51	0
unknown	SSC00003	-9	-9	1	-9	1.03	1.03	16945	11840	3495	159218	118.544	0.11	-139.27	135.12	1718.49	-1230.76	1844.47	0
unknown	SSC00004	-9	-9	1	-9	1.02	1.02	19159	11841	5762	157004	116.745	0.12	-160.11	-1011.12	2903.40	-320.78	3020.00	0
unknown	SSC00005	-9	-9	2	-9	1.99	0.00	17241	11827	3969	158922	118.284	0.11	1159.92	164.66	1301.97	-948.43	2678.17	0
unknown	SSC00006	-9	-9	1	-9	1.03	1.03	16743	11836	3273	159420	118.710	0.11	465.95	-1278.65	629.43	-1894.03	1565.34	0
unknown	SSC00011	-9	-9	1	-9	1.02	1.01	16761	11841	3292	159402	118.689	0.11	-647.15	1066.58	2230.69	-614.67	2766.30	0
"""

colors = sns.color_palette('Set2', 4)

for i, y in enumerate(sorted(df.sex.unique())):
    sub = df.sex == y
    ax.scatter(df['CNX'][sub], df['CNY'][sub], facecolors=colors[i],
            edgecolors=(0.8, 0.8, 0.8),
            label='inferred CN for chrX: %d' % y)


ax.set_xlabel("X copy number")
ax.set_ylabel("Y copy number")
ax.set_xticks([0, 1, 2])
ax.set_yticks([0, 1, 2])
legend = ax.legend(frameon=True, loc="upper left")
legend.get_frame().set_facecolor('#fafafa')

sns.despine()

plt.show()
