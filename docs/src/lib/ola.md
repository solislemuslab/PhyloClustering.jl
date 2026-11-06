# OLA Encoding of Phylogenetic Trees

OLA (ordered leaf attachment) encodes a rooted binary tree on `n` labelled leaves as a
vector of `n-1` integers: entry `a[i]` is the node that leaf `i+1` attaches next to when
the tree is built by adding leaves in order. Encoding and decoding both run in `O(n)`.

The valid codes are exactly those with `-i < a[i] ≤ i` and `a[i] ≠ 0` for every `i`.
Every such vector decodes to a tree, so the encoding is a bijection between rooted
binary trees and this set, and validity can be checked without decoding.

OLA needs a rooted tree, so [`ola_encode`](@ref) ignores the root of the input and
re-roots the tree on the pendant edge of an outgroup taxon. The code depends on this
choice and on the leaf ordering. Only the topology is encoded.

```@docs
ola_encode
ola_decode
ola_isvalid
ola_random_code
```
