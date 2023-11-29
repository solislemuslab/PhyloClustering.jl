# Hierarchical Clustering

[Hierarchical clustering](https://en.wikipedia.org/wiki/Hierarchical_clustering)
algorithms build a dendrogram of nested clusters by repeatedly merging
or splitting clusters.

The `hclust` function implements several classical algorithms for hierarchical
clustering (the algorithm to use is defined by the `linkage` parameter):

```@docs
hc_label
```

Single-linkage clustering using distance matrix:
```@example
using Clustering
D = rand(1000, 1000);
D += D'; # symmetric distance matrix (optional)
result = hclust(D, linkage=:single)
```
