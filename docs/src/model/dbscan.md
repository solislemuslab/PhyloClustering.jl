# DBSCAN

[Density-based Spatial Clustering of Applications with Noise
(DBSCAN)](http://en.wikipedia.org/wiki/DBSCAN) is a data clustering
algorithm that finds clusters through density-based expansion of seed
points. The algorithm was proposed in:

> Martin Ester, Hans-peter Kriegel, Jörg S, and Xiaowei Xu *A
> density-based algorithm for discovering clusters in large spatial
> databases with noise.* 1996.

## Interface

The implementation of *DBSCAN* algorithm provided by [Clustering](https://juliastats.org/Clustering.jl/stable/dbscan.html)
supports the two ways of specifying clustering data:

```@docs
dbscan_label
```
