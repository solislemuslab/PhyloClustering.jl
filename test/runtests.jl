using PhyloClustering
using Test

tests = [
    "test_bipartition.jl",
    "test_helper.jl",
    "test_models.jl"
]

anyerrors = false

for t in tests
    global anyerrors
    try
        @testset "$t" begin
            include(t)
        end
        println("\033[1m\033[32mPASSED\033[0m: $t")
    catch
        anyerrors = true
        println("\033[1m\033[31mFAILED\033[0m: $t")
    end
end
println("-------------------------------------")

if anyerrors
    throw("Tests failed")
else
    println("\033[1m\033[32mTests passed")
end
