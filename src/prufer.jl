# Prüfer sequences — labeled trees ↔ integer sequences of length n−2
#
# Note: "labeled" here means **all vertices** carry labels 1..n, not just the leaves.
# Phylogenetic trees name only their leaves, so this does not transfer directly: a
# labeling for the internal vertices has to be supplied first, which is what
# `net_to_adj` does — see reading/Prufer/notes.md.
#
# Key rule: vertex v appears in the sequence deg(v) − 1 times.
#   Appearing 0 times ⟺ v is a leaf. This single rule is both what encode produces
#   and what decode relies on.
#
# Complexity: the method itself is O(n) in both directions. The implementations here
#   are the naive O(n²) ones — each step scans linearly for the smallest leaf. A heap
#   gives O(n log n); a forward-only pointer gives O(n). See the TODO at the end.

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

Algorithm: give every vertex a degree — one, plus the number of times it appears in the
sequence. Then read the sequence left to right: the lowest-numbered vertex of degree one
is the current leaf, so join it to the entry being read and take one off both degrees.
When the sequence runs out, exactly two vertices of degree one are left; join them.

A tree on n vertices has n−1 edges but the sequence holds only n−2 numbers: that last
edge is free, since once the sequence is exhausted there is only one way to join the two
vertices that still need an edge.

Every sequence over 1..n of length n−2 maps to a tree (the map is surjective), so this
also serves as a **uniform random tree generator**: `prufer_decode(rand(1:n, n-2), n)`.
"""
function prufer_decode(seq::Vector{Int}, n::Int)
    length(seq) == n - 2 || throw(ArgumentError("sequence length must be n−2"))
    all(x -> 1 <= x <= n, seq) || throw(ArgumentError("sequence entries must be in 1..n"))
    n >= 2 || throw(ArgumentError("requires at least 2 vertices"))
    n == 2 && return [(1, 2)]

    # deg[v] = how many edges v still needs. By "v appears deg(v) − 1 times" this is one
    # more than the appearances v has left in the sequence, so deg[v] == 1 marks a leaf of
    # whatever is left. Attaching a vertex takes it to 0, which keeps it out of the scan
    # without a separate "used" flag.
    deg = ones(Int, n)
    for x in seq
        deg[x] += 1
    end

    edges = Vector{Tuple{Int,Int}}()
    sizehint!(edges, n - 1)

    for x in seq
        v = findfirst(==(1), deg)    # lowest-numbered leaf; this scan is the O(n²)
        v === nothing && error("no leaf found: this should not happen")

        push!(edges, minmax(v, x))
        deg[v] -= 1                  # v is attached and out of the running
        deg[x] -= 1                  # x used up one of the edges it still needed
    end

    # The sequence is exhausted; exactly two vertices still need an edge.
    rest = findall(==(1), deg)
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


# TODO Linear time. Both directions can reach O(n) with a forward-only pointer: after
#      deleting the smallest leaf v, its neighbour p may become a leaf — if p < v it is
#      the next one to go, so handle it in place and leave the pointer alone; otherwise
#      let the pointer reach it. The pointer only ever moves forward, so the scanning is
#      O(n) amortised over the whole run. This is a classical result, not new work; a
#      heap is the easier O(n log n) middle ground.
#
# TODO The image on leaf-labeled trees. `net_to_adj` supplies a canonical labeling, so
#      one tree gives one sequence and `prufer_encode(::HybridNetwork)` is well defined.
#      What is still open is the other direction: the map is injective but far from
#      surjective — at n = 6, 105 unrooted binary trees against 10^8 sequences — and no
#      simple description of the image (the validity conditions) is known. See the last
#      section of reading/Prufer/notes.md.
