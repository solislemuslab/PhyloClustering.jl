# Prüfer Sequence Encoding of Phylogenetic Trees

A Prüfer sequence encodes a tree on `n` labelled vertices as a sequence of `n-2`
integers. On fully labelled trees the encoding is a bijection, which is the standard
proof of Cayley's formula.

Phylogenetic trees are leaf-labelled, so the internal nodes have to be numbered before
the algorithm can run. [`net_to_adj`](@ref) does that in a way that depends only on the
tree itself, and not on how its Newick string happened to be written.

On leaf-labelled trees the encoding is injective but not surjective: most sequences of
the right length decode to a graph that is not a valid phylogenetic tree, so
[`prufer_decode`](@ref) is not the inverse of `prufer_encode(::HybridNetwork)`.

```@docs
net_to_adj
prufer_encode
prufer_decode
```
