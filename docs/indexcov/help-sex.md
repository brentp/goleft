Inferred Sex
============

Since we have scaled coverage to have a median of 1, chromosomes with a value of ~0.5 have a single
copy of that chromosome while those with a value of 1 have 2 copies (assuming diploid organism).

Using that information, we can infer the copy-number state of the sex chromosomes and report.

For example:
![Sex Help](https://cloud.githubusercontent.com/assets/1739/22120288/6b304990-de3c-11e6-811e-fb01b0a83a67.png "sex help")

In these plots **each point represents a sample**.

Here, we can see inferred males in blue and inferred females in pink. (Color choices aside...)

We can also see a single outlier sample that has only a single X chromosome and no Y chromosomes. In this case
we know that there is a problem with this sample where much of the data from the Y chromosome was lost.

*When viewed in context, this plot is interactive*; the name of the sample will appear in a tool-tip when the 
user hovers over a point.

Atypical Ploidy
---------------

This plot will also enable discover of samples with unusual sex chromosome copy-numbers. For example, an **XXY**
sample would appear at *2* along the X-axis and *1* along the Y-axis (there are no points there in this cohort).
An **XYY** sample would appear above the males in this plot at position *2* on the Y-axis. The occurence of This
in the population is as high as 1 in 400 so it should not be too unexpected to find such samples in large cohorts.

