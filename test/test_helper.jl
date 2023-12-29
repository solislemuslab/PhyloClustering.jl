include("../src/helper.jl")
@testset "helper.jl" begin
    @testset "standardize trees" begin
        trees = readMultiTopology("file/4-taxon-tree.trees")
        tree = split_weight(trees, 4)
        N = size(tree)[1]
        B = size(tree)[2]
        @test !(mean(tree[1,:]) ≈ 0)
        tree = standardize_tree(tree)
        @test size(tree) == (B,N)
        @test [std(tree[i,:]) for i in 1:7] ≈ [1,1,1,1,1,0,1] atol=0.1
        @test [mean(tree[i,:]) for i in 1:7] ≈ [0,0,0,0,0,0,0] atol=0.1

        trees = readMultiTopology("file/8-taxon-tree.trees")
        tree = split_weight(trees, 8)
        N = size(tree)[1]
        B = size(tree)[2]
        n = num_bipartitions(8)
        @test !(mean(tree[1,:]) ≈ 0)
        tree = standardize_tree(tree)
        @test size(tree) == (B,N)
        @test [std(tree[i,:]) for i in 1:8] ≈ repeat([1], outer = 8) atol=0.1
        @test [mean(tree[i,:]) for i in 1:n] ≈ repeat([0], outer = n) atol=0.1

        trees = readMultiTopology("file/16-taxon-tree.trees")
        tree = split_weight(trees, 16)
        N = size(tree)[1]
        B = size(tree)[2]
        n = num_bipartitions(16)
        @test !(mean(tree[1,:]) ≈ 0)
        tree = standardize_tree(tree)
        @test size(tree) == (B,N)
        @test [std(tree[i,:]) for i in 1:16] ≈ repeat([1], outer = 16) atol=0.1
        @test [mean(tree[i,:]) for i in 1:n] ≈ repeat([0], outer = n) atol=0.1
    end

    @testset "visualization" begin
        trees = readMultiTopology("file/4-taxon-tree.trees")
        tree = split_weight(trees, 4)
        label = [1,2,1,2]
        @test_logs plot_clusters(tree, label)
    end

    @testset "Euclidean distance matrix" begin
        trees = readMultiTopology("file/8-taxon-tree.trees")
        tree = split_weight(trees, 8)
        tree = standardize_tree(tree)
        matrix = distance(tree)
        @test size(matrix)[1] == size(matrix)[2]
        @test [matrix[i,i] for i in 1:size(matrix)[1]] == repeat([0], outer = size(matrix)[1])
        @test matrix ≈ [
            0.0 5.096169006818051 5.40722180738128 5.625211680819353; 
            5.096169006818051 0.0 4.37996068613932 5.896362019739273; 
            5.40722180738128 4.37996068613932 0.0 6.260740103673938; 
            5.625211680819353 5.896362019739273 6.260740103673938 0.0] atol=0.001

        trees = readMultiTopology("file/16-taxon-tree.trees")
        tree = split_weight(trees, 16)
        tree = standardize_tree(tree)
        matrix = distance(tree)
        @test size(matrix)[1] == size(matrix)[2]
        @test [matrix[i,i] for i in 1:size(matrix)[1]] == repeat([0], outer = size(matrix)[1])
        @test matrix ≈ [
            0.0 11.3979509761689 11.549754911257478 11.358576754022506; 
            11.3979509761689 0.0 11.21122183718288 11.930804141875372; 
            11.549754911257478 11.21122183718288 0.0 11.473318029384393; 
            11.358576754022506 11.930804141875372 11.473318029384393 0.0] atol=0.001
    end
end