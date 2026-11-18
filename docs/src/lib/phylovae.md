# PhyloVAE Encoding of Phylogenetic Trees

PhyloVAE encodes an unrooted binary tree on `n` labelled leaves as a vector of `n-3`
integers: entry `k` names the node that leaf `k+3` attaches next to when the tree is built
by adding leaves in order. Encoding and decoding both run in `O(n)`.

Nodes are named rather than edges, which is what keeps both directions linear: leaves take
`1:n` in order of addition, internal nodes take `(n+1):(2n-3)` in order of creation, and a
name is never reused or changed. The node joining the first three leaves is never named.

The valid codes are exactly those whose entry `k` lies in `1:(k+2)` or `(n+1):(n+k-1)`, so
validity can be checked without decoding. Unlike the single interval of
[`ola_isvalid`](@ref), each entry's set is two intervals: the gap `(k+3):n` holds the leaves
that have not been added yet. A coordinate-wise mean of valid codes can therefore land in a
gap, which is worth knowing before averaging codes.

This is section 3.2 of the paper. The variational autoencoder built on top of this encoding
is not implemented here.

Only the topology is encoded, and the code depends on the leaf ordering. No outgroup is
needed: the encoding is defined on unrooted trees, so [`phylovae_encode`](@ref) and
[`phylovae_decode`](@ref) are mutually inverse, where [`ola_encode`](@ref) first re-roots and
so is only a retraction.

```@docs
phylovae_encode
phylovae_decode
phylovae_isvalid
phylovae_random_code
```
