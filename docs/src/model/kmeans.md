# K-means

[K-means](http://en.wikipedia.org/wiki/K_means) is a classical method for
clustering or vector quantization. It produces a fixed number of clusters,
each associated with a *center* (also known as a *prototype*), and each data
point is assigned to a cluster with the nearest center.

From a mathematical standpoint, K-means is a coordinate descent
algorithm that solves the following optimization problem:
```math
\text{minimize} \ \sum_{i=1}^n \| \mathbf{x}_i - \boldsymbol{\mu}_{z_i} \|^2 \ \text{w.r.t.} \ (\boldsymbol{\mu}, z)
```
Here, ``\boldsymbol{\mu}_k`` is the center of the ``k``-th cluster, and
``z_i`` is an index of the cluster for ``i``-th point ``\mathbf{x}_i``.

```@docs
kmeans_label
```

## Examples

```@example
```