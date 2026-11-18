using PhyloNetworks: readnewick, readmultinewick, tiplabels
using PhyloClustering: _pv_build, _pv_code, _pv_postorder

# all valid codes of length n−3: entry k ranges over 1:(k+2) and (n+1):(n+k-1)
pv_codes(n) = Iterators.product(ntuple(k -> vcat(1:(k + 2), (n + 1):(n + k - 1)), n - 3)...)

# set of nontrivial splits: equal iff same unrooted topology
function pv_splits(s)
    _, c1, c2, c3, lab, root, N = _pv_build(s)
    below = Dict{Int,Set{Int}}()
    out = Set{Set{Int}}()
    for v in _pv_postorder(c1, c2, c3, root)
        if c1[v] == 0
            below[v] = Set([lab[v]])
        else
            below[v] = union(below[c1[v]], below[c2[v]])
            c3[v] != 0 && union!(below[v], below[c3[v]])
        end
        v == root && continue
        side = 1 in below[v] ? setdiff(Set(1:N), below[v]) : below[v]
        2 <= length(side) <= N - 2 && push!(out, side)
    end
    return out
end

doublefactorial(m) = prod(1:2:m)

@testset "phylovae.jl" begin

    # the paper's Figure 1 tree ↔ (3, 7, 5)
    @testset "worked example, n = 6" begin
        s = [3, 7, 5]
        tree = phylovae_decode(s)
        @test sort(tiplabels(tree)) == ["1", "2", "3", "4", "5", "6"]
        @test pv_splits(s) == Set([Set([3, 4]), Set([5, 6]), Set([3, 4, 5, 6])])
        @test phylovae_encode(tree) == s
        @test phylovae_encode(readnewick("((1,2),(5,6),(3,4));")) == s
    end

    # unlike the OLA box, each entry's valid set is two intervals with a gap
    @testset "the valid set has a gap" begin
        @test phylovae_isvalid([3, 7, 5])
        @test phylovae_isvalid([1, 1, 1])
        @test !phylovae_isvalid([4])            # n = 4, entry 1 must be in 1:3
        @test !phylovae_isvalid([3, 5, 5])      # n = 6, entry 2 must be in 1:4 or 7:7
        @test !phylovae_isvalid([3, 6, 5])      # 6 is a leaf not yet added at that step
        @test !phylovae_isvalid([3, 7, 9])      # n = 6, entry 3 must be in 1:5 or 7:8
        @test !phylovae_isvalid([0])
        @test !phylovae_isvalid(Int[])

        for n in 4:8
            @test length(collect(pv_codes(n))) == doublefactorial(2n - 5)
            @test all(phylovae_isvalid(collect(s)) for s in pv_codes(n))
        end

        # the gap at step k is exactly the leaves not yet added
        for n in 6:8, k in 1:(n - 3)
            gap = setdiff(1:(n + k - 1), vcat(1:(k + 2), (n + 1):(n + k - 1)))
            @test gap == collect((k + 3):n)
        end
    end

    @testset "round trip, code to tree and back" begin
        for n in 4:8
            @test all(_pv_code(_pv_build(collect(s))...) == collect(s) for s in pv_codes(n))
        end
    end

    @testset "exhaustive, (2n−5)!! distinct trees" begin
        for n in 4:7
            trees = Set{Set{Set{Int}}}()
            for s in pv_codes(n)
                push!(trees, pv_splits(collect(s)))
            end
            @test length(trees) == doublefactorial(2n - 5)     # 3, 15, 105, 945
        end
    end

    @testset "indexed internal nodes are exactly n+1..2n−3" begin
        # r is never indexed, so there are n−3 of them, not n−2
        for n in 4:7, s in pv_codes(n)
            _, c1, c2, c3, lab, root, N = _pv_build(collect(s))
            cf  = zeros(Int, length(c1))
            idx = Int[]
            for v in _pv_postorder(c1, c2, c3, root)
                if c1[v] == 0
                    cf[v] = lab[v]
                elseif v == root
                    cf[v] = min(cf[c1[v]], cf[c2[v]], cf[c3[v]])
                else
                    x, y = cf[c1[v]], cf[c2[v]]
                    cf[v] = min(x, y)
                    push!(idx, N + max(x, y) - 3)
                end
            end
            @test sort(idx) == collect((n + 1):(2n - 3))
        end
    end

    @testset "random codes are valid and decode" begin
        for n in 4:12, _ in 1:20
            s = phylovae_random_code(n)
            @test length(s) == n - 3
            @test phylovae_isvalid(s)
            @test length(tiplabels(phylovae_decode(s))) == n
        end
        @test_throws ArgumentError phylovae_random_code(3)
    end

    @testset "phylovae_encode on HybridNetwork" begin
        # the same unrooted tree written four ways (same fixtures as net_to_adj)
        same = ["(4,(1,(2,3)));",
                "(3,(2,(1,4)));",
                "(1,4,(2,3));",
                "(2,3,(4,1));"]
        codes = [phylovae_encode(readnewick(s)) for s in same]

        @testset "independent of how the Newick is written" begin
            @test all(==(codes[1]), codes)
            @test codes[1] == [1]
        end

        @testset "independent of branch lengths" begin
            @test phylovae_encode(readnewick("(4:4.2,(1:2.5,(2:2.1,3:2.1):0.4):1.8);")) == codes[1]
        end

        @testset "leaves numbered by sorted name" begin
            @test phylovae_encode(readnewick("(D,(A,(B,C)));")) == codes[1]
        end

        @testset "a different topology gives a different code" begin
            @test phylovae_encode(readnewick("(1,2,(3,4));")) != codes[1]
        end

        # no re-rooting and no outgroup, so unlike ola_encode this is a true inverse
        @testset "encode and decode are mutually inverse" begin
            for s in vcat(same, ["((1,2),(5,6),(3,4));", "((1,3),2,(4,5));",
                                 "(5,(1,(2,(3,4))));"])
                net = readnewick(s)
                a = phylovae_encode(net)
                @test phylovae_encode(phylovae_decode(a)) == a
            end
            for n in 4:7, s in pv_codes(n)
                s = collect(s)
                @test phylovae_encode(phylovae_decode(s)) == s
            end
        end

        @testset "six taxa, from the test fixtures" begin
            # the file also holds a tree with a polytomy, which must be rejected
            for t in readmultinewick("file/6-taxon-tree.trees")
                if all(length(nd.edge) <= 3 for nd in t.node)
                    s = phylovae_encode(t)
                    @test length(s) == 3
                    @test phylovae_isvalid(s)
                    @test phylovae_encode(phylovae_decode(s)) == s
                else
                    @test_throws ArgumentError phylovae_encode(t)
                end
            end
        end

        @testset "argument is left alone" begin
            net = readnewick("(4,(1,(2,3)));")
            before = length(net.node)
            phylovae_encode(net)
            @test length(net.node) == before
        end
    end

    @testset "input validation" begin
        @test_throws ArgumentError phylovae_encode(readnewick("(1,2,3,4,5);"))  # not binary
        @test_throws ArgumentError phylovae_encode(readnewick("(1,2,3);"))      # n < 4
        @test_throws ArgumentError phylovae_decode([4])                          # outside the set
        @test_throws ArgumentError phylovae_decode(Int[])                        # n < 4
        @test_throws ArgumentError phylovae_decode([3]; taxa = ["A", "B"])       # wrong count
        @test_throws ArgumentError phylovae_decode([3]; taxa = ["A", "A", "B", "C"])
    end

end
