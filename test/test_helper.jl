include("../src/helper.jl")
@testset "helper.jl" begin
    @testset "standardize trees" begin
        trees = readmultinewick("file/4-taxon-tree.trees")
        tree = split_weight(trees, 4)
        N = size(tree)[1]
        B = size(tree)[2]
        @test !(mean(tree[1, :]) ≈ 0)
        tree = standardize_tree(tree)
        @test size(tree) == (B, N)
        @test [std(tree[i, :]) for i in 1:7] ≈ [1, 1, 1, 1, 1, 0, 1] atol = 0.1
        @test [mean(tree[i, :]) for i in 1:7] ≈ [0, 0, 0, 0, 0, 0, 0] atol = 0.1

        trees = readmultinewick("file/8-taxon-tree.trees")
        tree = split_weight(trees, 8)
        N = size(tree)[1]
        B = size(tree)[2]
        n = num_bipartitions(8)
        @test !(mean(tree[1, :]) ≈ 0)
        tree = standardize_tree(tree)
        @test size(tree) == (B, N)
        @test [std(tree[i, :]) for i in 1:8] ≈ repeat([1], outer=8) atol = 0.1
        @test [mean(tree[i, :]) for i in 1:n] ≈ repeat([0], outer=n) atol = 0.1

        trees = readmultinewick("file/16-taxon-tree.trees")
        tree = split_weight(trees, 16)
        N = size(tree)[1]
        B = size(tree)[2]
        n = num_bipartitions(16)
        @test !(mean(tree[1, :]) ≈ 0)
        tree = standardize_tree(tree)
        @test size(tree) == (B, N)
        @test [std(tree[i, :]) for i in 1:16] ≈ repeat([1], outer=16) atol = 0.1
        @test [mean(tree[i, :]) for i in 1:n] ≈ repeat([0], outer=n) atol = 0.1
    end

    @testset "visualization" begin
        trees = readmultinewick("file/4-taxon-tree.trees")
        tree = split_weight(trees, 4)
        label = [1, 2, 1, 2]
    end

    @testset "Euclidean distance matrix" begin
        trees = readmultinewick("file/8-taxon-tree.trees")
        tree = split_weight(trees, 8)
        tree = standardize_tree(tree)
        matrix = distance(tree)
        @test size(matrix)[1] == size(matrix)[2]
        @test [matrix[i, i] for i in 1:size(matrix)[1]] == repeat([0], outer=size(matrix)[1])
        @test matrix ≈ [
            0.0 5.096169006818051 5.40722180738128 5.625211680819353;
            5.096169006818051 0.0 4.37996068613932 5.896362019739273;
            5.40722180738128 4.37996068613932 0.0 6.260740103673938;
            5.625211680819353 5.896362019739273 6.260740103673938 0.0] atol = 0.001

        trees = readmultinewick("file/16-taxon-tree.trees")
        tree = split_weight(trees, 16)
        tree = standardize_tree(tree)
        matrix = distance(tree)
        @test size(matrix)[1] == size(matrix)[2]
        @test [matrix[i, i] for i in 1:size(matrix)[1]] == repeat([0], outer=size(matrix)[1])
        @test matrix ≈ [
            0.0 11.568782625290618 11.718373756037687 11.566452183632407;
            11.568782625290618 0.0 11.21122183718288 11.917964061705202;
            11.718373756037687 11.21122183718288 0.0 11.459965379864949;
            11.566452183632407 11.917964061705202 11.459965379864949 0.0] atol = 0.001
    end
end