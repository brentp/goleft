Depth
=====

Depth plots are very intuitive. The x-value indicates the position on the chromosome and
the y-value indicates the relative coverage.

For example:

![depth Help](https://cloud.githubusercontent.com/assets/1739/22121191/4505f720-de40-11e6-8f92-8e386298c2cd.png "depth help")

In this plot **each line represents a sample**.

Here, we can see that, as expected, most samples are centered around 1. We can see the centromere. Notably, we also see
a single sample with extremely low coverage. We can tell this is not due to a deletion because it is non-zero and it's 
not 0.5 (which would indicate a heterozygous deletion). For this case, we know that this occured due to a processing
'glitch' that resulted in a loss of most of the data for this region in that sample.

We can see a corresponding drop in the coverage plot below.

The depth plot in the `index.html` file is a static image. But clicking it will take the user to an interactive
version of that plot for more in-depth (so to speak) exploration.

Coverage
========

Coverage plots indicate the proportion of bins (on the y) covered to at least a relative depth (on the x).
We expect the large drop around 1 since much of each chromosome will have an estimated, scaled coverage of 1.

For example:

![coverage Help](https://cloud.githubusercontent.com/assets/1739/22121534/bf3e7cdc-de41-11e6-9486-1e81f8fbb2d8.png "coverage help")

In this plot **each line represents a sample**.

Where we can see that most samples drop around a coverage value (x-axis) of 1. There is a single
sample that drops much sooner. This sample is also seen in the depth plot above.

In general, changes as severe as this indicate a problem with the data rather than a biological effect such as a
deletion. Though a very large (20MB+) deletion would show up in this type of plot.

The coverage plot in the `index.html` file is a static image. But clicking it will take the user to an interactive
version of that plot for more in-depth exploration.
