# Shared tree representations and tools. Five encodings are built on top of this.
#
# Keep both adjacency list and edge list representations                                                                            
#   (1) adjacency list for encoding. (degrees/neighbors) adj::Vector{Vector{Int}} 
#   (2) edge list for decoding. (push) edges::Vector{Tuple{Int,Int}}
# Conversions are provided in both directions. When this is hooked up to PhyloClustering.jl,
# its tree type is a different one — this is where the adapter layer belongs.

"""
    edges_to_adj(edges, n) -> Vector{Vector{Int}}

Convert an edge list to an adjacency list. Vertex labels are assumed to be 1..n.
"""
function edges_to_adj(edges::Vector{Tuple{Int,Int}}, n::Int)
    adj = [Int[] for _ in 1:n]
    for (u, v) in edges
        push!(adj[u], v)
        push!(adj[v], u)
    end
    return adj
end


"""
    adj_to_edges(adj) -> Vector{Tuple{Int,Int}}

Convert an adjacency list to an edge list. Each edge appears once, as (u, v) with u < v,
and the result is sorted so that two trees can be compared directly.
"""
function adj_to_edges(adj::Vector{Vector{Int}})
    edges = Tuple{Int,Int}[]
    for u in eachindex(adj), v in adj[u]
        u < v && push!(edges, (u, v))
    end
    return sort(edges)
end


"""
    degrees(adj) -> Vector{Int}
"""
degrees(adj::Vector{Vector{Int}}) = length.(adj) #broadcast to get degrees of each vertex


"""
    leaves(adj) -> Vector{Int}

Vertices of degree 1, in ascending order.
"""
leaves(adj::Vector{Vector{Int}}) = [v for v in eachindex(adj) if length(adj[v]) == 1]


"""
    is_tree(edges, n) -> Bool

n vertices, n−1 edges, connected ⟹ a tree. Connectivity is checked with union-find.
"""
function is_tree(edges::Vector{Tuple{Int,Int}}, n::Int)
    length(edges) == n - 1 || return false
    parent = collect(1:n)
    find(x) = parent[x] == x ? x : (parent[x] = find(parent[x]))
    for (u, v) in edges
        ru, rv = find(u), find(v)
        ru == rv && return false      # cycle
        parent[ru] = rv
    end
    return true
end


"""
    random_tree(n) -> Vector{Tuple{Int,Int}}

Draw a uniformly random tree on n labeled vertices; returns an edge list.

Method: generate a random Prüfer sequence of length n−2 and decode it. Since Prüfer is a
bijection, a uniform sequence gives a uniform tree — this is the bijection put to work.
"""
function random_tree(n::Int)
    n >= 2 || throw(ArgumentError("Requires at least 2 vertices"))
    n == 2 && return [(1, 2)]
    return prufer_decode(rand(1:n, n - 2), n)
end

"""
    net_to_adj(net::HybridNetwork) -> Vector{Vector{Int}}

Convert a phylogenetic tree into a fully labeled adjacency list.

Phylogenetic trees are leaf-labeled: only the tips carry names. Prüfer coding needs
every vertex to carry a label, so the internal nodes have to be numbered too, and that
numbering has to depend on the tree alone — not on how the Newick string happened to
be written.

# Numbering
 - Leaves take `1:n`, by position in `sort(tiplabels(net))`. This is the rule
   [`get_bipartition`](@ref) uses, so every encoding in this package agrees on leaf
   numbers; [`num_to_name`](@ref) recovers the original names.
 - Internal nodes take `n+1` and up, in the order a traversal first reaches them. The
   traversal starts at leaf 1 and, at each node, visits neighbours in increasing order
   of the smallest leaf number below them.

Leaf 1 is fixed by the taxon names, whereas the root of a `HybridNetwork` read from
Newick is an artifact of where the string happened to start — so starting the
traversal at leaf 1 is what makes the numbering canonical. Degree-2 nodes, the usual
form that artifact takes, are removed first, so that the same unrooted tree always
gives the same adjacency list however it was written.

The argument is not modified. Reticulations are not supported.

# Example
```jldoctest
julia> net = readnewick("(4:4.249,(1:2.457,(2:2.064,3:2.064):0.393):1.793);");

julia> length(net_to_adj(net))   # 4 leaves + 2 internal nodes
6
```
"""
function net_to_adj(net::HybridNetwork)
    net.numhybrids == 0 || throw(ArgumentError(
        "net_to_adj expects a tree, but this network has $(net.numhybrids) reticulation(s)"))

    work = deepcopy(net)          # the caller's network is left alone
    removedegree2nodes!(work)     # drop the root artifact, if there is one

    taxa = sort(tiplabels(work))
    n = length(taxa)
    n >= 3 || throw(ArgumentError("requires at least 3 taxa"))
    leafnum = Dict(name => i for (i, name) in enumerate(taxa))

    nodes = work.node
    V = length(nodes)
    pos = IdDict(nd => i for (i, nd) in enumerate(nodes))

    # adjacency over working indices, before any renumbering
    work_adj = [Int[] for _ in 1:V]
    for e in work.edge
        a, b = pos[e.node[1]], pos[e.node[2]]
        push!(work_adj[a], b)
        push!(work_adj[b], a)
    end

    # leaves get their final numbers immediately; internal nodes stay at 0 for now
    label = zeros(Int, V)
    for (i, nd) in enumerate(nodes)
        nd.leaf || continue
        haskey(leafnum, nd.name) ||
            throw(ArgumentError("leaf with no usable name: \"$(nd.name)\""))
        label[i] = leafnum[nd.name]
    end

    start = findfirst(==(1), label)
    start === nothing && error("leaf 1 not found: this should not happen")

    # root the tree at leaf 1 and record the order nodes are reached
    parent = zeros(Int, V)
    parent[start] = -1
    order = Int[]
    sizehint!(order, V)
    stack = [start]
    while !isempty(stack)
        u = pop!(stack)
        push!(order, u)
        for v in work_adj[u]
            v == parent[u] && continue
            parent[v] = u
            push!(stack, v)
        end
    end
    length(order) == V || throw(ArgumentError("the input tree is not connected"))

    # smallest leaf number in each subtree. Walking `order` backwards visits every
    # child before its parent, so one pass is enough.
    minleaf = fill(typemax(Int), V)
    for u in Iterators.reverse(order)
        label[u] == 0 || (minleaf[u] = min(minleaf[u], label[u]))
        p = parent[u]
        p > 0 && (minleaf[p] = min(minleaf[p], minleaf[u]))
    end

    # canonical traversal: at each node take the subtree holding the smallest leaf
    # first. Internal nodes are numbered in the order they come up.
    next = n + 1
    stack = [start]
    while !isempty(stack)
        u = pop!(stack)
        if label[u] == 0
            label[u] = next
            next += 1
        end
        kids = [v for v in work_adj[u] if v != parent[u]]
        sort!(kids, by = v -> minleaf[v], rev = true)  # reversed: smallest ends on top
        append!(stack, kids)
    end
    next == V + 1 || error("internal numbering missed a node: this should not happen")

    adj = [Int[] for _ in 1:V]
    for u in 1:V, v in work_adj[u]
        push!(adj[label[u]], label[v])
    end
    foreach(sort!, adj)
    return adj
end
