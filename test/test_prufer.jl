using PhyloNetworks: readnewick, readmultinewick
using PhyloClustering: edges_to_adj, adj_to_edges, degrees, leaves, is_tree, random_tree

@testset "prufer.jl" begin

    @testset "worked example, n = 6" begin
        edges = [(1, 4), (2, 4), (3, 4), (4, 5), (5, 6)]
        adj = edges_to_adj(edges, 6)
        @test prufer_encode(adj) == [4, 4, 4, 5]
        @test prufer_decode([4, 4, 4, 5], 6) == sort(edges)
    end

    # Surjectivity: every sequence decodes to a tree, and encoding it returns the
    # sequence unchanged.
    @testset "round trip, sequence to tree and back" begin
        for n in 4:9, _ in 1:50
            seq = rand(1:n, n - 2)
            adj = edges_to_adj(prufer_decode(seq, n), n)
            @test prufer_encode(adj) == seq
        end
    end

    # Injectivity: distinct trees cannot share a sequence.
    @testset "round trip, tree to sequence and back" begin
        for n in 4:9, _ in 1:50
            edges = random_tree(n)
            adj = edges_to_adj(edges, n)
            @test prufer_decode(prufer_encode(adj), n) == edges
        end
    end

    # Exhaustive for small n: decoding every sequence recovers Cayley's formula,
    # n^(n-2) distinct trees. Any deviation in the implementation shows up as a
    # duplicate or a missing tree, so this tests the bijection itself rather than
    # a handful of fixed outputs.
    @testset "exhaustive, Cayley's formula" begin
        for n in 4:6
            trees = Set{Vector{Tuple{Int,Int}}}()
            for seq in Iterators.product(ntuple(_ -> 1:n, n - 2)...)
                edges = prufer_decode(collect(seq), n)
                @test is_tree(edges, n)
                push!(trees, edges)
            end
            @test length(trees) == n^(n - 2)    # 16, 125, 1296
        end
    end

    @testset "occurrences equal degree minus one" begin
        for n in 4:9, _ in 1:20
            adj = edges_to_adj(random_tree(n), n)
            seq = prufer_encode(adj)
            for v in 1:n
                @test count(==(v), seq) == length(adj[v]) - 1
            end
            # the corollary: labels that never appear are exactly the leaves
            @test sort([v for v in 1:n if count(==(v), seq) == 0]) == leaves(adj)
        end
    end

    @testset "star tree gives a constant sequence" begin
        for n in 4:8, c in 1:n
            edges = sort([minmax(c, u) for u in 1:n if u != c])
            @test prufer_encode(edges_to_adj(edges, n)) == fill(c, n - 2)
        end
    end

    @testset "label n survives to the end" begin
        for n in 4:9, _ in 1:50
            adj = edges_to_adj(random_tree(n), n)
            # n is removed only if it is the smallest leaf at some step, which cannot
            # happen while any other leaf is still standing
            @test count(==(n), prufer_encode(adj)) == length(adj[n]) - 1
        end
    end

    @testset "net_to_adj" begin
        # the same unrooted tree, written four ways: a degree-2 root in the first
        # two, a trifurcation in the others, and different starting taxa throughout
        same = ["(4,(1,(2,3)));",
                "(3,(2,(1,4)));",
                "(1,4,(2,3));",
                "(2,3,(4,1));"]
        adjs = [net_to_adj(readnewick(s)) for s in same]

        @testset "shape" begin
            adj = adjs[1]
            @test length(adj) == 6                      # 4 leaves + 2 internal nodes
            @test all(length(adj[v]) == 1 for v in 1:4) # leaves take 1:n
            @test all(length(adj[v]) == 3 for v in 5:6) # internal nodes are binary
            @test is_tree(adj_to_edges(adj), 6)
        end

        # The point of the canonical numbering: the Newick string's root is an
        # artifact of where the string starts, and must not reach the output.
        @testset "independent of how the Newick is written" begin
            for a in adjs[2:end]
                @test a == adjs[1]
            end
            seqs = [prufer_encode(readnewick(s)) for s in same]
            @test all(==(seqs[1]), seqs)
            @test seqs[1] == [5, 6, 6, 5]
        end

        @testset "independent of branch lengths" begin
            @test net_to_adj(readnewick("(4:4.2,(1:2.5,(2:2.1,3:2.1):0.4):1.8);")) == adjs[1]
        end

        # Leaves are numbered by position in sort(tiplabels), so names that sort the
        # same way describe the same numbered tree.
        @testset "leaves numbered by sorted name" begin
            @test net_to_adj(readnewick("(D,(A,(B,C)));")) == adjs[1]
            @test prufer_encode(readnewick("(D,(A,(B,C)));")) == [5, 6, 6, 5]
        end

        @testset "a different topology gives a different sequence" begin
            @test prufer_encode(readnewick("(1,2,(3,4));")) != [5, 6, 6, 5]
        end

        @testset "six taxa, from the test fixtures" begin
            for t in readmultinewick("file/6-taxon-tree.trees")
                adj = net_to_adj(t)
                @test is_tree(adj_to_edges(adj), length(adj))
                @test all(length(adj[v]) == 1 for v in 1:6)
                @test length(prufer_encode(t)) == length(adj) - 2
            end
        end

        @testset "argument is left alone" begin
            net = readnewick("(4,(1,(2,3)));")
            before = length(net.node)
            net_to_adj(net)
            @test length(net.node) == before
        end
    end

    # TODO Small experiment: change one element of a sequence and measure how far the
    #      tree moves under the Robinson-Foulds distance. This is the evidence for the
    #      "does the distance mean anything geometrically" column of the comparison
    #      table.
end
