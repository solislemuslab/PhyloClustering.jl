# OLA (ordered leaf attachment): rooted binary trees on n leaves ↔ integer vectors of
# length n−1. Both directions O(n).
#
# Indexing: the paper is 0-based; here leaves are 1..n and internal nodes −1..−(n−1).
#   So a[1] = 1, the valid set is −i < a[i] ≤ i (a[i] ≠ 0), and a node's name is −(CS − 1).
#   Internal node −i is stored at array index n+i (see _ola_idx).
#
# Reference: Richman, Zhang, Matsen, "Vector encoding of phylogenetic trees by
# ordered leaf attachment".

using PhyloNetworks, Random

_ola_idx(x::Int, n::Int) = x > 0 ? x : n - x


"""
    ola_isvalid(a) -> Bool

Check whether `a` is a valid OLA code, i.e. `-i < a[i] ≤ i` and `a[i] ≠ 0` for every `i`.

# Example
```jldoctest
julia> ola_isvalid([1, -1, 2, -3])
true

julia> ola_isvalid([1, 3])
false
```
"""
function ola_isvalid(a::AbstractVector{Int})
    isempty(a) && return false
    for (i, v) in enumerate(a)
        (v != 0 && -i < v <= i) || return false
    end
    return true
end


"""
    ola_random_code(n; rng = Random.default_rng()) -> Vector{Int}

Draw a uniformly random OLA code on `n` leaves. Decoding it with [`ola_decode`](@ref)
gives a uniformly random rooted binary tree.
"""
function ola_random_code(n::Int; rng::Random.AbstractRNG = Random.default_rng())
    n >= 3 || throw(ArgumentError("requires at least 3 leaves, got $n"))
    a = Vector{Int}(undef, n - 1)
    for i in 1:(n - 1)
        v = rand(rng, 1:(2i - 1))           # 2i−1 choices: leaves 1..i, internal −1..−(i−1)
        a[i] = v <= i ? v : i - v
    end
    return a
end


# Trees are kept as flat arrays: par, c1, c2 (children), lab (leaf label, 0 if internal).
# A node is a leaf iff c1[v] == 0.

# Decode: at step i, insert node −i on the parent edge of a[i] and hang leaf i+1 below it.
function _ola_build(a::AbstractVector{Int})
    ola_isvalid(a) || throw(ArgumentError("not a valid OLA code: $a"))
    n = length(a) + 1
    N = 2n - 1

    par = zeros(Int, N)
    c1  = zeros(Int, N)
    c2  = zeros(Int, N)
    lab = zeros(Int, N)
    for v in 1:n
        lab[v] = v
    end

    root = n + 1                        # start from −1(leaf 1, leaf 2), since a[1] = 1
    c1[root], c2[root] = 1, 2
    par[1] = par[2] = root

    for i in 2:(n - 1)
        v = _ola_idx(a[i], n)
        w = n + i                       # new node −i
        l = i + 1                       # new leaf
        p = par[v]
        if p == 0
            root = w
        elseif c1[p] == v
            c1[p] = w
        else
            c2[p] = w
        end
        par[w] = p
        c1[w], par[v] = v, w
        c2[w], par[l] = l, w
    end

    return par, c1, c2, lab, root, n
end


# Iterative, to avoid stack overflow on deep (caterpillar) trees.
function _ola_postorder(c1::Vector{Int}, c2::Vector{Int}, root::Int)
    order = Int[]
    sizehint!(order, length(c1))
    stack = [root]
    while !isempty(stack)
        u = pop!(stack)
        push!(order, u)
        if c1[u] != 0
            push!(stack, c1[u])
            push!(stack, c2[u])
        end
    end
    return reverse!(order)
end


# Encode, in two passes:
#   1. postorder: clade-founder CF = min leaf below; name = −(max of children's CF − 1)
#   2. for i = n−1 down to 1: a[i] = label of leaf i+1's sister, then prune leaf i+1
function _ola_code(par::Vector{Int}, c1::Vector{Int}, c2::Vector{Int},
                   lab::Vector{Int}, root::Int, n::Int)
    N = length(c1)
    order = _ola_postorder(c1, c2, root)
    length(order) == 2n - 1 || throw(ArgumentError(
        "the tree has $(length(order)) nodes reachable from the root, expected $(2n - 1)"))

    cf   = zeros(Int, N)
    name = zeros(Int, N)
    for v in order
        if c1[v] == 0
            cf[v] = lab[v]
        else
            x, y = cf[c1[v]], cf[c2[v]]
            cf[v]   = min(x, y)
            name[v] = -(max(x, y) - 1)
        end
    end

    leafidx = zeros(Int, n)
    for v in order
        c1[v] == 0 && (leafidx[lab[v]] = v)
    end

    P, C1, C2 = copy(par), copy(c1), copy(c2)
    a = Vector{Int}(undef, n - 1)
    for i in (n - 1):-1:1
        v = leafidx[i + 1]
        p = P[v]
        s = C1[p] == v ? C2[p] : C1[p]          # sister
        a[i] = c1[s] == 0 ? lab[s] : name[s]     # names never change, so no recomputation
        g = P[p]
        if g != 0
            C1[g] == p ? (C1[g] = s) : (C2[g] = s)
        end
        P[s] = g
    end
    return a
end


"""
    ola_decode(a; taxa = string.(1:length(a)+1)) -> HybridNetwork

Decode an OLA code into a rooted binary tree. Leaf `j` is named `taxa[j]`; by default
`"j"`, matching the leaf numbering of [`get_bipartition`](@ref). Branch lengths are left
unset.

# Example
```jldoctest
julia> tree = ola_decode([1, -1, 2, -3]);

julia> length(tiplabels(tree))
5

julia> ola_encode(tree; outgroup = "3") == [1, -1, 2, -3]
true

julia> sort(tiplabels(ola_decode([1, -1, 2, -3]; taxa = ["A","B","C","D","E"])))
5-element Vector{String}:
 "A"
 "B"
 "C"
 "D"
 "E"
```
"""
function ola_decode(a::AbstractVector{Int}; taxa::AbstractVector{<:AbstractString} = string.(1:(length(a) + 1)))
    length(a) >= 2 || throw(ArgumentError("requires a code of length at least 2 (n ≥ 3)"))
    _, c1, c2, lab, root, n = _ola_build(a)
    length(taxa) == n ||
        throw(ArgumentError("got $(length(taxa)) taxon names for $n leaves"))
    allunique(taxa) || throw(ArgumentError("taxon names must be unique"))
    return readnewick(_ola_newick(c1, c2, lab, root, taxa))
end


# Iterative Newick writer. Sentinels on the stack: 0 = comma, −1 = closing paren.
function _ola_newick(c1::Vector{Int}, c2::Vector{Int}, lab::Vector{Int}, root::Int,
                     taxa::AbstractVector{<:AbstractString})
    io = IOBuffer()
    stack = [root]
    while !isempty(stack)
        u = pop!(stack)
        if u == 0
            print(io, ',')
        elseif u == -1
            print(io, ')')
        elseif c1[u] == 0
            print(io, taxa[lab[u]])
        else
            print(io, '(')
            push!(stack, -1)
            push!(stack, c2[u])
            push!(stack, 0)
            push!(stack, c1[u])
        end
    end
    print(io, ';')
    return String(take!(io))
end


"""
    ola_encode(tree::HybridNetwork; outgroup = :first) -> Vector{Int}

Encode a binary phylogenetic tree as an OLA code of length `n-1`, where `n` is the
number of leaves. Only the topology is encoded.

Leaves are numbered by position in `sort(tiplabels(tree))`, as in
[`get_bipartition`](@ref). The root of the input is ignored: the tree is treated as
unrooted and re-rooted on the pendant edge of `outgroup` (`:first` for the first taxon,
or a taxon name). Different outgroups generally give different codes.

# Example
```jldoctest
julia> ola_encode(readnewick("(1,2,(3,4));"))
3-element Vector{Int64}:
 1
 2
 3

julia> ola_encode(readnewick("(1,2,(3,4));"); outgroup = "3")
3-element Vector{Int64}:
  1
 -1
 -1

julia> ola_encode(readnewick("((1,((2,4),5)),3);"))
4-element Vector{Int64}:
  1
  2
  2
 -3
```
"""
function ola_encode(tree::HybridNetwork; outgroup::Union{Symbol,AbstractString} = :first)
    par, c1, c2, lab, root, n = _ola_root(tree, outgroup)
    return _ola_code(par, c1, c2, lab, root, n)
end


# HybridNetwork -> flat arrays, rooted on the pendant edge of the outgroup leaf.
function _ola_root(net::HybridNetwork, outgroup::Union{Symbol,AbstractString})
    net.numhybrids == 0 || throw(ArgumentError(
        "ola_encode expects a tree, but this network has $(net.numhybrids) reticulation(s)"))

    work = deepcopy(net)
    removedegree2nodes!(work)         # drop the Newick root, so the input's rooting is ignored

    taxa = sort(tiplabels(work))
    n = length(taxa)
    n >= 3 || throw(ArgumentError("requires at least 3 taxa, got $n"))
    leafnum = Dict(name => i for (i, name) in enumerate(taxa))

    nodes = work.node
    V = length(nodes)
    V == 2n - 2 || throw(ArgumentError(
        "expected an unrooted binary tree on $n leaves ($(2n-2) nodes), got $V nodes"))
    pos = IdDict(nd => i for (i, nd) in enumerate(nodes))

    adj = [Int[] for _ in 1:V]
    for e in work.edge
        x, y = pos[e.node[1]], pos[e.node[2]]
        push!(adj[x], y)
        push!(adj[y], x)
    end

    lab0 = zeros(Int, V)
    for (i, nd) in enumerate(nodes)
        nd.leaf || continue
        haskey(leafnum, nd.name) ||
            throw(ArgumentError("leaf with no usable name: \"$(nd.name)\""))
        lab0[i] = leafnum[nd.name]
    end
    for i in 1:V
        d = length(adj[i])
        if lab0[i] == 0
            d == 3 || throw(ArgumentError(
                "not binary: an internal node has degree $d, expected 3"))
        else
            d == 1 || throw(ArgumentError(
                "leaf \"$(nodes[i].name)\" has degree $d, expected 1"))
        end
    end

    ogname = outgroup === :first ? taxa[1] : String(outgroup)
    haskey(leafnum, ogname) ||
        throw(ArgumentError("outgroup \"$ogname\" is not a taxon of this tree"))
    ℓ = findfirst(==(leafnum[ogname]), lab0)

    N = V + 1
    R = N                             # new root, on the pendant edge of ℓ
    par = zeros(Int, N)
    c1  = zeros(Int, N)
    c2  = zeros(Int, N)
    lab = zeros(Int, N)
    lab[1:V] .= lab0

    nb = only(adj[ℓ])
    c1[R], c2[R] = ℓ, nb
    par[ℓ] = R
    par[nb] = R

    # `from` tracks the adjacency neighbour we came from; par[nb] is R, which is not in adj
    from = zeros(Int, V)
    from[nb] = ℓ
    stack = [nb]
    while !isempty(stack)
        u = pop!(stack)
        kids = [v for v in adj[u] if v != from[u]]
        isempty(kids) && continue
        length(kids) == 2 || throw(ArgumentError(
            "not binary once rooted: a node has $(length(kids)) children"))
        c1[u], c2[u] = kids[1], kids[2]
        for v in kids
            par[v] = u
            from[v] = u
            push!(stack, v)
        end
    end

    return par, c1, c2, lab, R, n
end
