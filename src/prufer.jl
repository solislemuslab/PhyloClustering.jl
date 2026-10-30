# Prüfer sequences — labeled trees ↔ integer sequences of length n−2
#
# Note: "labeled" here means **all vertices** carry labels 1..n, not just the leaves.
# Phylogenetic trees are leaf-labeled with anonymous internal nodes, so this does not
# transfer directly — see reading/Prufer/notes.md.
#
# Key rule: vertex v appears in the sequence deg(v) − 1 times.
#   Appearing 0 times ⟺ v is a leaf. This single rule is both what encode produces
#   and what decode relies on.
#
# Complexity: the naive implementation is O(n²) — each step scans linearly for the
#   smallest leaf / candidate. A heap brings this down to O(n log n); see the TODO
#   at the end of the file.

"""
    prufer_encode(adj) -> Vector{Int}

Encode a tree on n labeled vertices as a sequence of length n−2.

Algorithm: repeat n−2 times — find the smallest-labeled leaf v, append v's only
neighbour to the sequence, and remove v.

# Example
```jldoctest
julia> adj = PhyloClustering.edges_to_adj([(1,4),(2,4),(3,4),(4,5),(5,6)], 6);

julia> prufer_encode(adj)
4-element Vector{Int64}:
 4
 4
 4
 5
```
"""
function prufer_encode(adj::Vector{Vector{Int}})
    n = length(adj)
    n >= 2 || throw(ArgumentError("requires at least 2 vertices"))
    n == 2 && return Int[]

    # Do not touch the caller's adj: degrees and the "removed" flags are kept separately.
    deg = length.(adj)
    removed = falses(n)
    seq = Vector{Int}(undef, n - 2)

    for k in 1:(n - 2)
        # Smallest leaf. Naive linear scan — this is where the O(n²) comes from.
        v = 0
        for u in 1:n
            if !removed[u] && deg[u] == 1
                v = u
                break                # scanning upward, so the first hit is the smallest
            end
        end
        v == 0 && error("no leaf found: input is not a tree")

        # v is a leaf, so it has exactly one surviving neighbour.
        nb = 0
        for u in adj[v]
            if !removed[u]
                nb = u
                break
            end
        end
        nb == 0 && error("leaf has no neighbour: input is not a tree")

        seq[k] = nb
        removed[v] = true
        deg[v] = 0
        deg[nb] -= 1                 # nb lost a neighbour, and may now be a leaf itself
    end

    return seq
end


"""
    prufer_decode(seq, n) -> Vector{Tuple{Int,Int}}

Decode a Prüfer sequence of length n−2 into a tree on n labeled vertices, returned as
an edge list (each edge as (u, v) with u < v, sorted globally).

Algorithm: maintain the candidate leaves — labels absent from the **remaining** sequence
and not yet used. At each step take the smallest candidate, join it to the first element
of the sequence, then consume that element; if the consumed label no longer appears in
the remaining sequence, it becomes a candidate in turn. Once the sequence is exhausted,
join the two labels that are left.

Every sequence over 1..n of length n−2 maps to a tree (the map is surjective), so this
also serves as a **uniform random tree generator**: `prufer_decode(rand(1:n, n-2), n)`.
"""
function prufer_decode(seq::Vector{Int}, n::Int)
    length(seq) == n - 2 || throw(ArgumentError("sequence length must be n−2"))
    all(x -> 1 <= x <= n, seq) || throw(ArgumentError("sequence entries must be in 1..n"))
    n >= 2 || throw(ArgumentError("requires at least 2 vertices"))
    n == 2 && return [(1, 2)]

    # count[v] = how many times v still appears in the remaining sequence. Initialised
    # from the whole sequence. By "appearances = deg − 1", count[v] == 0 means v is
    # currently a leaf.
    count = zeros(Int, n)
    for x in seq
        count[x] += 1
    end

    used = falses(n)                        # vertices already attached and removed from the tree
    edges = Vector{Tuple{Int,Int}}()
    sizehint!(edges, n - 1)

    for x in seq
        # Smallest candidate leaf: not yet used, and absent from the remaining sequence.
        v = 0
        for u in 1:n
            if !used[u] && count[u] == 0
                v = u
                break
            end
        end
        v == 0 && error("no candidate leaf found: this should not happen")

        push!(edges, minmax(v, x))
        used[v] = true
        count[x] -= 1                       # consume this appearance; at 0, x becomes a candidate
    end

    # The sequence is exhausted; exactly two vertices are left, one of which must be n.
    rest = [u for u in 1:n if !used[u]]
    length(rest) == 2 || error("not exactly 2 vertices left: this should not happen")
    push!(edges, (rest[1], rest[2]))

    return sort(edges)
end

"""
    prufer_encode(tree::HybridNetwork) -> Vector{Int}

Encode a phylogenetic tree, by way of [`net_to_adj`](@ref).

An unrooted binary tree on n leaves has 2n−2 vertices, so the sequence has length
2n−4. Unlike the labeled-tree case this map is injective but **not** surjective: most
sequences of that length decode to a graph that is not a valid phylogenetic tree,
since nothing in the sequence forces the internal vertices to degree 3 or the leaf
labels to degree 1. `prufer_decode` is therefore not its inverse. Characterising the
image — the validity conditions — is still open; see `reading/Prufer/notes.md`.
"""
prufer_encode(tree::HybridNetwork) = prufer_encode(net_to_adj(tree))


# TODO Heap optimisation: maintain the candidate leaves in a min-heap, bringing both
#      encode and decode to O(n log n). Get the naive version working first; record
#      O(n²) in the comparison table and note that it can be improved.
#
# TODO Apply this to leaf-labeled trees — the question the project actually needs
#      answered. Phylogenetic trees have unlabeled internal nodes, so a canonical
#      labeling is required first; that labeling is not unique, which means the same
#      tree can yield different sequences and the bijection breaks at that level.
#      See the last section of reading/Prufer/notes.md.
