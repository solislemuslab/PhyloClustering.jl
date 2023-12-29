# K-means

[K-means](http://en.wikipedia.org/wiki/K_means) is a classical method for
clustering or vector quantization. It produces a fixed number of clusters,
each associated with a *center* (also known as a *prototype*), and each data
point is assigned to a cluster with the nearest center.

The K-means clustering used in our package is Yinyang K-means which has less runtime and memory usage on large datasets. Traditional K-means calculates the distance from all data points to centroids for each iteration. Yinyang K-means uses triangular inequalities to construct and maintain upper and lower bounds on the distances of data points from the centroids, with global and local filtering to minimize unneeded calculations.

```@docs
kmeans_label
```

## [Example](@id usage)

```@example 1
using PhyloClustering, PhyloNetworks, StableRNGs
# use RNG for stable result
rng = StableRNG(1)

# read trees with 4-taxa in Newick format using PhyloNetworks
trees = readMultiTopology("../data/data.trees");

# convert trees to Bipartition foramt and embed them via split-weight embedding
trees = split_weight(trees, 4)
```

Standardize the data and input them into Yinyang K-means clustering.

```@example 1
tree = standardize_tree(trees);
label = kmeans_label(tree, 2, rng=rng)
```

We can visualize the result using build-in function [`plot_clusters`](@ref).

```@example 1
plot_clusters(trees', label)
```

**Reference**

The implementation of *Yinyang K-means* is provided by [`ParallelKMeans.jl`](https://github.com/PyDataBlog/ParallelKMeans.jl).

> Yufei Ding et al. 2015. *Yinyang K-Means: A Drop-In Replacement of the Classic K-Means with Consistent Speedup*. 
> Proceedings of the 32nd International Conference on Machine Learning, ICML 2015, Lille, France, 6-11 July 2015