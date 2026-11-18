# PhyloVAE encoding: unrooted binary trees on N leaves ↔ integer vectors of length N−3.
# Both directions O(N). This is section 3.2 of the paper only; the VAE is not implemented.
#
# Build τ₃ on leaves 1,2,3 around an unlabelled node r, then add leaf n+1 by splitting the
# parent edge of the node indexed s[n−2]. Leaves are indexed 1..N, internal nodes N+1..2N−3
# in creation order, and r is never indexed — hence N−3 indexed internal nodes, not N−2.
#
# Reference: Xie, Richman, Gao, Matsen IV, Zhang, "PhyloVAE: unsupervised learning of
# phylogenetic trees via variational autoencoders", ICLR 2025, section 3.2.

using PhyloNetworks, Random


"""
    phylovae_isvalid(s) -> Bool

Check whether `s` is a valid PhyloVAE code: entry `k` must index a node present in the tree
of `k+2` leaves, i.e. lie in `1:(k+2)` (a leaf) or `(N+1):(N+k-1)` (an internal node), where
`N = length(s) + 3`.

The valid set is a product of per-entry sets, so validity is checked without decoding. Unlike
[`ola_isvalid`](@ref), each set is two intervals rather than one: the gap `(k+3):N` holds the
leaves that have not been added yet.

# Example
```jldoctest
julia> phylovae_isvalid([3, 7, 5])
true

julia> phylovae_isvalid([3, 5, 5])
false
```
"""
function phylovae_isvalid(s::AbstractVector{Int})
    isempty(s) && return false
    N = length(s) + 3
    for (k, v) in enumerate(s)
        (1 <= v <= k + 2 || N + 1 <= v <= N + k - 1) || return false
    end
    return true
end


"""
    phylovae_random_code(n; rng = Random.default_rng()) -> Vector{Int}

Draw a uniformly random PhyloVAE code on `n` leaves. Decoding it with
[`phylovae_decode`](@ref) gives a uniformly random unrooted binary tree, since the code set
is in bijection with those trees and has the same size `(2n-5)!!`.
"""
function phylovae_random_code(n::Int; rng::Random.AbstractRNG = Random.default_rng())
    n >= 4 || throw(ArgumentError("requires at least 4 leaves, got $n"))
    s = Vector{Int}(undef, n - 3)
    for k in 1:(n - 3)
        v = rand(rng, 1:(2k + 1))        # 2k+1 = 2(k+2)−3 choices at this step
        s[k] = v <= k + 2 ? v : n + v - k - 2
    end
    return s
end


# Trees are kept as flat arrays: par, c1, c2, c3, lab (leaf label, 0 if internal).
# c3 is nonzero only at r, which has three children; a node is a leaf iff c1[v] == 0.

# Decode: at step k, insert node N+k on the parent edge of s[k] and hang leaf k+3 below it.
function _pv_build(s::AbstractVector{Int})
    phylovae_isvalid(s) || throw(ArgumentError("not a valid PhyloVAE code: $s"))
    N = length(s) + 3
    A = 2N - 2

    par = zeros(Int, A)
    c1  = zeros(Int, A)
    c2  = zeros(Int, A)
    c3  = zeros(Int, A)
    lab = zeros(Int, A)
    for v in 1:N
        lab[v] = v
    end

    r = A
    c1[r], c2[r], c3[r] = 1, 2, 3
    par[1] = par[2] = par[3] = r

    for k in 1:(N - 3)
        v = s[k]
        t = N + k                       # new internal node, index N+n−2 with n = k+2
        l = k + 3                       # new leaf
        p = par[v]
        if c1[p] == v
            c1[p] = t
        elseif c2[p] == v
            c2[p] = t
        else
            c3[p] = t
        end
        par[t] = p
        c1[t], par[v] = v, t
        c2[t], par[l] = l, t
    end

    return par, c1, c2, c3, lab, r, N
end


# Iterative, to avoid stack overflow on deep (caterpillar) trees.
function _pv_postorder(c1::Vector{Int}, c2::Vector{Int}, c3::Vector{Int}, root::Int)
    order = Int[]
    sizehint!(order, length(c1))
    stack = [root]
    while !isempty(stack)
        u = pop!(stack)
        push!(order, u)
        if c1[u] != 0
            push!(stack, c1[u])
            push!(stack, c2[u])
            c3[u] != 0 && push!(stack, c3[u])
        end
    end
    return reverse!(order)
end


# Encode, in two passes:
#   1. postorder: cf = min leaf below; an internal node ≠ r is indexed N + max(children cf) − 3,
#      because it was created when that leaf arrived and split the path below it
#   2. for n = N−1 down to 3: s[n−2] = index of leaf n+1's sister, then prune leaf n+1
function _pv_code(par::Vector{Int}, c1::Vector{Int}, c2::Vector{Int}, c3::Vector{Int},
                  lab::Vector{Int}, root::Int, N::Int)
    A = length(c1)
    order = _pv_postorder(c1, c2, c3, root)
    length(order) == 2N - 2 || throw(ArgumentError(
        "the tree has $(length(order)) nodes reachable from r, expected $(2N - 2)"))

    cf  = zeros(Int, A)
    idx = zeros(Int, A)
    for v in order
        if c1[v] == 0
            cf[v] = lab[v]
            idx[v] = lab[v]
        elseif v == root
            cf[v] = min(cf[c1[v]], cf[c2[v]], cf[c3[v]])
        else
            x, y = cf[c1[v]], cf[c2[v]]
            cf[v]  = min(x, y)
            idx[v] = N + max(x, y) - 3
        end
    end

    leafpos = zeros(Int, N)
    for v in order
        c1[v] == 0 && (leafpos[lab[v]] = v)
    end

    P, C1, C2, C3 = copy(par), copy(c1), copy(c2), copy(c3)
    s = Vector{Int}(undef, N - 3)
    for n in (N - 1):-1:3
        v = leafpos[n + 1]
        p = P[v]
        sis = C1[p] == v ? C2[p] : C1[p]        # p is never r, so it has exactly two children
        s[n - 2] = idx[sis]                      # indices never change, so no recomputation
        g = P[p]
        if C1[g] == p
            C1[g] = sis
        elseif C2[g] == p
            C2[g] = sis
        else
            C3[g] = sis
        end
        P[sis] = g
    end
    return s
end


"""
    phylovae_decode(s; taxa = string.(1:length(s)+3)) -> HybridNetwork

Decode a PhyloVAE code into an unrooted binary tree. Leaf `j` is named `taxa[j]`; by default
`"j"`, matching the leaf numbering of [`get_bipartition`](@ref). Branch lengths are left
unset.

# Example
```jldoctest
julia> tree = phylovae_decode([3, 7, 5]);

julia> length(tiplabels(tree))
6

julia> phylovae_encode(tree) == [3, 7, 5]
true

julia> sort(tiplabels(phylovae_decode([1, 2]; taxa = ["A","B","C","D","E"])))
5-element Vector{String}:
 "A"
 "B"
 "C"
 "D"
 "E"
```
"""
function phylovae_decode(s::AbstractVector{Int};
                         taxa::AbstractVector{<:AbstractString} = string.(1:(length(s) + 3)))
    length(s) >= 1 || throw(ArgumentError("requires a code of length at least 1 (n ≥ 4)"))
    _, c1, c2, c3, lab, r, N = _pv_build(s)
    length(taxa) == N ||
        throw(ArgumentError("got $(length(taxa)) taxon names for $N leaves"))
    allunique(taxa) || throw(ArgumentError("taxon names must be unique"))
    return readnewick(_pv_newick(c1, c2, c3, lab, r, taxa))
end


# Iterative Newick writer. Sentinels on the stack: 0 = comma, −1 = closing paren.
# r is written as a trifurcation, so the output is unrooted.
function _pv_newick(c1::Vector{Int}, c2::Vector{Int}, c3::Vector{Int}, lab::Vector{Int},
                    root::Int, taxa::AbstractVector{<:AbstractString})
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
            if c3[u] != 0
                push!(stack, c3[u])
                push!(stack, 0)
            end
            push!(stack, c2[u])
            push!(stack, 0)
            push!(stack, c1[u])
        end
    end
    print(io, ';')
    return String(take!(io))
end


"""
    phylovae_encode(tree::HybridNetwork) -> Vector{Int}

Encode a binary phylogenetic tree as a PhyloVAE code of length `n-3`, where `n` is the number
of leaves. Only the topology is encoded.

Leaves are numbered by position in `sort(tiplabels(tree))`, as in [`get_bipartition`](@ref).
The code depends on that ordering. The root of the input is ignored and no outgroup is
needed: the encoding is defined on unrooted trees, with the node joining the first three
leaves playing the role that a root plays in [`ola_encode`](@ref).

# Example
```jldoctest
julia> phylovae_encode(readnewick("((1,2),(5,6),(3,4));"))
3-element Vector{Int64}:
 3
 7
 5

julia> phylovae_encode(readnewick("(1,2,(3,4));"))
1-element Vector{Int64}:
 3

julia> phylovae_encode(readnewick("((1,3),2,(4,5));"))
2-element Vector{Int64}:
 2
 4
```
"""
function phylovae_encode(tree::HybridNetwork)
    par, c1, c2, c3, lab, r, N = _pv_root(tree)
    return _pv_code(par, c1, c2, c3, lab, r, N)
end


# HybridNetwork -> flat arrays, rooted at r = the node where the paths to leaves 1, 2 and 3 meet.
function _pv_root(net::HybridNetwork)
    net.numhybrids == 0 || throw(ArgumentError(
        "phylovae_encode expects a tree, but this network has $(net.numhybrids) reticulation(s)"))

    work = deepcopy(net)
    removedegree2nodes!(work)         # drop the Newick root, so the input's rooting is ignored

    taxa = sort(tiplabels(work))
    N = length(taxa)
    N >= 4 || throw(ArgumentError("requires at least 4 taxa, got $N"))
    leafnum = Dict(name => i for (i, name) in enumerate(taxa))

    nodes = work.node
    V = length(nodes)
    V == 2N - 2 || throw(ArgumentError(
        "expected an unrooted binary tree on $N leaves ($(2N-2) nodes), got $V nodes"))
    pos = IdDict(nd => i for (i, nd) in enumerate(nodes))

    adj = [Int[] for _ in 1:V]
    for e in work.edge
        x, y = pos[e.node[1]], pos[e.node[2]]
        push!(adj[x], y)
        push!(adj[y], x)
    end

    lab = zeros(Int, V)
    for (i, nd) in enumerate(nodes)
        nd.leaf || continue
        haskey(leafnum, nd.name) ||
            throw(ArgumentError("leaf with no usable name: \"$(nd.name)\""))
        lab[i] = leafnum[nd.name]
    end
    for i in 1:V
        d = length(adj[i])
        if lab[i] == 0
            d == 3 || throw(ArgumentError(
                "not binary: an internal node has degree $d, expected 3"))
        else
            d == 1 || throw(ArgumentError(
                "leaf \"$(nodes[i].name)\" has degree $d, expected 1"))
        end
    end

    leafpos = zeros(Int, N)
    for i in 1:V
        lab[i] != 0 && (leafpos[lab[i]] = i)
    end
    r = _pv_median(adj, leafpos[1], leafpos[2], leafpos[3])

    par = zeros(Int, V)
    c1  = zeros(Int, V)
    c2  = zeros(Int, V)
    c3  = zeros(Int, V)
    stack = [r]
    while !isempty(stack)
        u = pop!(stack)
        kids = [v for v in adj[u] if v != par[u]]
        isempty(kids) && continue
        nk = u == r ? 3 : 2
        length(kids) == nk || throw(ArgumentError(
            "not binary once rooted at r: a node has $(length(kids)) children, expected $nk"))
        c1[u], c2[u] = kids[1], kids[2]
        u == r && (c3[u] = kids[3])
        for v in kids
            par[v] = u
            push!(stack, v)
        end
    end

    return par, c1, c2, c3, lab, r, N
end


# The node where the paths to leaves a, b, c meet: root at a, then take the LCA of b and c.
function _pv_median(adj::Vector{Vector{Int}}, a::Int, b::Int, c::Int)
    V = length(adj)
    up = zeros(Int, V)
    seen = falses(V)
    seen[a] = true
    stack = [a]
    while !isempty(stack)
        u = pop!(stack)
        for v in adj[u]
            if !seen[v]
                seen[v] = true
                up[v] = u
                push!(stack, v)
            end
        end
    end

    onpath = falses(V)
    u = b
    while u != 0
        onpath[u] = true
        u == a && break
        u = up[u]
    end
    u = c
    while !onpath[u]
        u = up[u]
    end
    return u
end
