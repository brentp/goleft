Principal Components Analysis
=============================

In `indexcov` we use PCA as a means of dimensionality reduction. That is, we send a matrix of $n samples by $m bins in
the genome and (with PCA), we project that to $n samples by 5 components.

This enables discovery of major batch effects among samples.
**NOTE** that `indexcov` excludes the sex chromosomes from the PCA, otherwise the 1st or 2nd principal Components
would separate the sexes because of coverage difference on the sex chromosomes. Because `indexcov` already examines
sex in the [sex plot](https://github.com/brentp/goleft/blob/master/docs/indexcov/help-sex.md) it is excluded from the
PCA so that we can more easily discover other batch effects.

We plot the *1st vs 2nd* principal components:

![pca1](https://cloud.githubusercontent.com/assets/1739/22120760/54e89e38-de3e-11e6-8bc0-c0ee22f23f99.png "pca1")


and the *1st vs 3rd* principal components:

![pca2](https://cloud.githubusercontent.com/assets/1739/22120759/54d8aa5a-de3e-11e6-8b46-f0e6c1d94080.png "pca2")

In this case, we don't see any major clusterings or outliers, but we may want to hover over the samples at the
edges and check their sample Ids in the other `indexcov` plots.
