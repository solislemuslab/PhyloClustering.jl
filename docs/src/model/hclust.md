# Hierarchical Clustering

[Hierarchical clustering](https://en.wikipedia.org/wiki/Hierarchical_clustering) is a hierarchical method that creates a clustering tree called a dendrogram. Each leaf on the tree is a data point, and branches represent the joining of clusters. We can 'cut' the tree at different heights to get a different number of clusters. This method does not require the number of clusters to be specified in advance and can be either agglomerative (bottom-up) or divisive (top-down).

```@docs
hc_label
```
## Example
```@example 1
using PhyloClustering, PhyloNetworks

# read trees with 4-taxa in Newick format using PhyloNetworks
trees = readMultiTopology("../data/data.trees");

# convert trees to Bipartition foramt and embed them via split-weight embedding
trees = split_weight(trees, 4);
```

Standardize the data, calculate the distance matrix, and input them into Yinyang K-means clustering.

```@example 1
tree = standardize_tree(trees);
matrix = distance(tree);
label = hc_label(matrix, 2)
```

We can visualize the result using build-in function [`plot_clusters`](@ref).

```@example 1
plot_clusters(trees', label)
```

**Reference**

The implementation of *hierarchical clustering* is provided by [`Clustering.jl`](https://github.com/JuliaStats/Clustering.jl).

> Joe H. Ward Jr. (1963) *Hierarchical Grouping to Optimize an Objective Function*, 
> Journal of the American Statistical Association, 58:301, 236-244, DOI: 10.1080/01621459.1963.10500845
