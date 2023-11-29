include("../src/models.jl")
using StableRNGs
using Suppressor

rng = StableRNG(376)

trees = readMultiTopology("file/8-taxon-tree.trees")
tree = print_bipartition(trees, 8)
tree_8 = standardize_tree(tree)

trees = readMultiTopology("file/16-taxon-tree.trees")
tree = print_bipartition(trees, 16)
tree_16 = standardize_tree(tree);

@testset "models.jl" begin
    @testset "Yinyang K-means clustering" begin
        @test kmeans_label(tree_8, 1) == [1, 1, 1, 1]

        Random.seed!(rng, 376)
        @test kmeans_label(tree_8, 2, rng=rng) == [1, 1, 2, 1]
        @test kmeans_label(tree_8, 2, init="rand", rng=rng) == [1, 2, 2, 1]
        @test kmeans_label(tree_16, 3, rng=rng) == [3, 1, 2, 3]
    end

    @testset "hierarchical clustering" begin
        matrix = distance(tree_16)
        @test hc_label(matrix, 1) == [1, 1, 1, 1]
        @test hc_label(matrix, 2) == [1, 2, 2, 1]

        matrix = distance(tree_8)
        @test hc_label(matrix, 2) == [1, 1, 1, 2]
        @test hc_label(matrix, 3, linkage=:average) == [1, 2, 2, 3]
        @test hc_label(matrix, 2, linkage=:single) == [1, 1, 1, 2]
    end

    @testset "gaussian mixture model" begin
        @suppress begin
            matrix = distance(tree_8)
            @test gmm_label(matrix, 2, method=:split) == ([2, 2, 2, 1], [2, 2, 2, 1])
            @test gmm_label(tree_8, 2, method=:split) == ([1, 1, 1, 1], [1, 1, 1, 1])
            @test gmm_label(tree_16, 2, method=:split) == ([1, 1, 1, 1], [1, 1, 1, 1])
        end
    end

    @testset "Density-based Spatial Clustering of Applications with Noise" begin
        Random.seed!(376)
        X = hcat(
            randn(2, 5) .+ [0., 5.],
            randn(2, 5) .+ [-5., 0.],
            randn(2, 5) .+ [5., 0.])
        @test dbscan_label(X, 3) == [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3]

        Random.seed!(376)
        X = randn(2, 15)
        @test dbscan_label(X,0.9,min_cluster_size=3) == [1, 0, 1, 1, 1, 2, 1, 2, 1, 0, 1, 1, 2, 2, 1]
        @test dbscan_label(X,0.8,min_neighbors=2) == [1, 0, 1, 1, 1, 2, 1, 3, 1, 0, 1, 1, 3, 2, 1]
    end
end