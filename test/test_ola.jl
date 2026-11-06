using PhyloNetworks: readnewick, readmultinewick, tiplabels
using PhyloClustering: _ola_build, _ola_code, _ola_postorder

# all valid codes of length n−1
ola_box(n) = Iterators.product(ntuple(i -> [v for v in -(i - 1):i if v != 0], n - 1)...)

# set of clades: equal iff same rooted topology
function ola_clades(a)
    _, c1, c2, lab, root, _ = _ola_build(a)
    seen = Dict{Int,Set{Int}}()
    out = Set{Set{Int}}()
    for v in _ola_postorder(c1, c2, root)
        s = c1[v] == 0 ? Set([lab[v]]) : union(seen[c1[v]], seen[c2[v]])
        seen[v] = s
        push!(out, s)
    end
    return out
end

doublefactorial(m) = prod(1:2:m)

@testset "ola.jl" begin

    # paper's (0,−1,1,−3) ↔ ((0,((1,3),4)),2), shifted to 1-based
    @testset "worked example, n = 5" begin
        a = [1, -1, 2, -3]
        tree = ola_decode(a)
        @test sort(tiplabels(tree)) == ["1", "2", "3", "4", "5"]
        @test ola_clades(a) == Set([Set([1]), Set([2]), Set([3]), Set([4]), Set([5]),
                                    Set([2, 4]), Set([2, 4, 5]), Set([1, 2, 4, 5]),
                                    Set([1, 2, 3, 4, 5])])
        @test ola_encode(tree; outgroup = "3") == a
    end

    @testset "the valid set is a box" begin
        @test ola_isvalid([1, -1, 2, -3])
        @test ola_isvalid([1, 1, 1, 1])
        @test !ola_isvalid([1, 3])          # a[2] must satisfy −2 < a[2] ≤ 2
        @test !ola_isvalid([1, -2])
        @test !ola_isvalid([1, 0])          # 0 is not a node label
        @test !ola_isvalid([2])             # a[1] can only be 1
        @test !ola_isvalid(Int[])

        for n in 3:7
            @test length(collect(ola_box(n))) == doublefactorial(2n - 3)
            @test all(ola_isvalid(collect(a)) for a in ola_box(n))
        end
    end

    @testset "round trip, code to tree and back" begin
        for n in 3:7
            @test all(_ola_code(_ola_build(collect(a))...) == collect(a) for a in ola_box(n))
        end
    end

    @testset "exhaustive, (2n−3)!! distinct trees" begin
        for n in 3:6
            trees = Set{Set{Set{Int}}}()
            for a in ola_box(n)
                push!(trees, ola_clades(collect(a)))
            end
            @test length(trees) == doublefactorial(2n - 3)     # 3, 15, 105, 945
        end
    end

    @testset "a[1] is always 1" begin
        for n in 3:7
            @test all(collect(a)[1] == 1 for a in ola_box(n))
        end
        for s in ["(1,2,(3,4));", "((1,((2,4),5)),3);", "(5,(1,(2,(3,4))));"]
            @test ola_encode(readnewick(s))[1] == 1
        end
    end

    @testset "internal names are exactly −1..−(n−1)" begin
        # Proposition 6
        for n in 3:6, a in ola_box(n)
            _, c1, c2, lab, root, _ = _ola_build(collect(a))
            names = Int[]
            cf = zeros(Int, length(c1))
            for v in _ola_postorder(c1, c2, root)
                if c1[v] == 0
                    cf[v] = lab[v]
                else
                    x, y = cf[c1[v]], cf[c2[v]]
                    cf[v] = min(x, y)
                    push!(names, -(max(x, y) - 1))
                end
            end
            @test sort(names) == collect(-(n - 1):-1)
        end
    end

    # Proposition 8
    @testset "last-leaf deletion is stable" begin
        for n in 4:6, a in ola_box(n)
            a = collect(a)
            par, c1, c2, lab, root, _ = _ola_build(a)
            v = n
            p = par[v]
            s = c1[p] == v ? c2[p] : c1[p]
            g = par[p]
            if g == 0
                root = s
            else
                c1[g] == p ? (c1[g] = s) : (c2[g] = s)
            end
            par[s] = g
            @test _ola_code(par, c1, c2, lab, root, n - 1) == a[1:end - 1]
        end
    end

    @testset "random codes are valid and decode" begin
        for n in 3:12, _ in 1:20
            a = ola_random_code(n)
            @test length(a) == n - 1
            @test ola_isvalid(a)
            @test length(tiplabels(ola_decode(a))) == n
        end
        @test_throws ArgumentError ola_random_code(2)
    end

    @testset "ola_encode on HybridNetwork" begin
        # the same unrooted tree written four ways (same fixtures as net_to_adj)
        same = ["(4,(1,(2,3)));",
                "(3,(2,(1,4)));",
                "(1,4,(2,3));",
                "(2,3,(4,1));"]
        codes = [ola_encode(readnewick(s)) for s in same]

        @testset "independent of how the Newick is written" begin
            @test all(==(codes[1]), codes)
            @test codes[1] == [1, 2, -2]
        end

        @testset "independent of branch lengths" begin
            @test ola_encode(readnewick("(4:4.2,(1:2.5,(2:2.1,3:2.1):0.4):1.8);")) == codes[1]
        end

        @testset "leaves numbered by sorted name" begin
            @test ola_encode(readnewick("(D,(A,(B,C)));")) == codes[1]
        end

        @testset "a different topology gives a different code" begin
            @test ola_encode(readnewick("(1,2,(3,4));")) != codes[1]
        end

        @testset "the outgroup changes the code" begin
            net = readnewick("((1,((2,4),5)),3);")
            @test ola_encode(net) == [1, 2, 2, -3]
            @test ola_encode(net; outgroup = "1") == [1, 2, 2, -3]
            @test ola_encode(net; outgroup = "3") == [1, -1, 2, -3]
            @test ola_encode(net; outgroup = "2") != ola_encode(net)
        end

        # encode re-roots, so it is only a retraction: encode ∘ decode ∘ encode = encode
        @testset "encode is a retraction" begin
            for s in vcat(same, ["((1,((2,4),5)),3);", "(5,(1,(2,(3,4))));"])
                a = ola_encode(readnewick(s))
                @test ola_encode(ola_decode(a)) == a
            end
        end

        @testset "six taxa, from the test fixtures" begin
            # the file also holds a tree with a polytomy, which must be rejected
            for t in readmultinewick("file/6-taxon-tree.trees")
                if all(length(nd.edge) <= 3 for nd in t.node)
                    a = ola_encode(t)
                    @test length(a) == 5
                    @test ola_isvalid(a)
                    @test ola_encode(ola_decode(a)) == a
                else
                    @test_throws ArgumentError ola_encode(t)
                end
            end
        end

        @testset "argument is left alone" begin
            net = readnewick("(4,(1,(2,3)));")
            before = length(net.node)
            ola_encode(net)
            @test length(net.node) == before
        end
    end

    @testset "input validation" begin
        @test_throws ArgumentError ola_encode(readnewick("(1,2,3,4);"))       # not binary
        @test_throws ArgumentError ola_encode(readnewick("(1,2,(3,4));"); outgroup = "99")
        @test_throws ArgumentError ola_decode([1, 3])                          # outside the box
        @test_throws ArgumentError ola_decode([1])                             # n < 3
        @test_throws ArgumentError ola_decode([1, -1]; taxa = ["A", "B"])      # wrong count
        @test_throws ArgumentError ola_decode([1, -1]; taxa = ["A", "A", "B"]) # duplicates
    end

end
