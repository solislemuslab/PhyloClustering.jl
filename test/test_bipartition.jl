include("../src/bipartition.jl")
using Suppressor # for capture printed result

@testset "bipartition.jl" begin
    @testset "show bipartitions of certain number of taxa or between certain indices" begin
        @test num_bipartitions(4) == 7
        @test num_bipartitions(8) == 127
        @test num_bipartitions(16) == 32767
        @test (@capture_out show_bipartitions(4)) == "idx\tpartition\n0\t1 | 2 3 4 \n1\t2 | 1 3 4 \n2\t3 | 1 2 4 \n3\t4 | 1 2 3 \n4\t1 2 | 3 4 \n5\t1 3 | 2 4 \n6\t1 4 | 2 3 \n"
        @test (@capture_out show_bipartitions(7, start = 20, stop = 26)) == "idx\tpartition\n20\t3 6 | 1 2 4 5 7 \n21\t3 7 | 1 2 4 5 6 \n22\t4 5 | 1 2 3 6 7 \n23\t4 6 | 1 2 3 5 7 \n24\t4 7 | 1 2 3 5 6 \n25\t5 6 | 1 2 3 4 7 \n26\t5 7 | 1 2 3 4 6 \n"
        @test (@capture_out show_bipartitions(3, start = 5, stop = 3)) == "idx\tpartition\n"
        @test (@capture_out show_bipartition(10, 277)) == "idx\tpartition\n277\t 2  3  8  9 |  1  4  5  6  7 10 \n"
    end

    @testset "get bipartitions of phylogenetic trees in Newick format " begin
        trees = readMultiTopology("file/4-taxon-tree.trees")
        @test get_bipartition(trees[2],4) == [3 => 4.04, 2 => 2.035, 0 => 1.084, 1 => 1.084, 4 => 0.95, 3 => 2.006]
        @test num_to_name(trees[1]) == Dict{Int64, Any}(4 => "4", 2 => "2", 3 => "3", 1 => "1")
        @test (@capture_out get_bipartition(trees[1],6)) == "check the parameter n"
        @test split_weight(trees,4) == [
            2.073 2.073 2.492 7.764 0.419 0.0 0.0; 
            1.084 1.084 2.035 6.045999999999999 0.95 0.0 0.0; 
            1.234 1.234 2.1 6.221 0.865 0.0 0.0; 
            2.457 2.064 2.064 6.042 0.0 0.0 0.393]

        trees = readMultiTopology("file/test.trees")
        @test get_bipartition(trees[4],4) == [1 => 3.555816, 2 => 2.399569, 0 => 0.629763, 3 => 0.629763, 6 => 1.769806, 1 => 1.156247]
        @test num_to_name(trees[1]) == Dict{Int64, Any}(4 => "P2", 2 => "O", 3 => "P1", 1 => "HYB")
        @test split_weight(trees,4) == [
            1.593392 6.139258 1.207652 1.207652 0.38574 0.0 0.0; 
            0.852455 2.093894 2.093894 0.852455 0.0 0.0 1.779641; 
            1.534422 2.006191 2.287441 1.534422 0.0 0.0 0.471769; 
            0.629763 4.712063000000001 2.399569 0.629763 0.0 0.0 1.769806; 
            1.524985 4.543632 1.524985 2.4089 0.0 0.883915 0.0]

        trees = readMultiTopology("file/6-taxon-tree.trees")
        @test get_bipartition(trees[3],6) == [5 => 0.08797164, 4 => 0.07148751, 2 => 0.04760093, 3 => 0.05347179, 0 => 0.0342379, 1 => 0.02864877, 6 => 0.005502266, 22 => 0.01711993, 6 => 0.004750963]
        @test (@capture_out get_bipartition(trees[2],11)) == "check the parameter n"
        @test split_weight(trees, 6) == [
            0.0256169 0.01151126 0.03644864 0.03587151 0.07877169 0.06931551 0.03743503 0.0 0.0 0.0 0.0 0.003646167 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.01414154 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
            0.01862177 0.04619028 0.02224014 0.02355061 0.084584 0.1004835 0.03322833 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.01557455 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
            0.0342379 0.02864877 0.04760093 0.05347179 0.07148751 0.08797164 0.010253228999999999 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.01711993 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 
            0.02527073 0.01873581 0.02161674 0.0488751 0.06761764 0.09904998 0.030438164 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.03232537 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    end
end