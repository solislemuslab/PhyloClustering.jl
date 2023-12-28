using PhyloNetworks, Combinatorics, Formatting

"""
        num_bipartitions(n::Int64)

    Function to get number of bipartitions.

    # Arguments
        - `n`: the number of taxa.
    # Output:
        The number of bipartition corresponding to the number of taxa.
"""
function num_bipartitions(n::Int64)
    return 2^(n-1) - 1
end;


"""
        show_bipartitions(n::Int64; start::Int64 = 0, stop::Int64=-1)

    Function to show the bipartitions between certain indices of the number of taxa.

    # Arguments
        `n`: the number of taxa
        `start`(defaults to `0`): the starting index.
        `stop`(defaults to `-1`): the ending index.
    # Output
        The bipartiions of trees with n taxa between start index and stop index.
    # Example
    show_bipartitions(4,stop = 4):
        idx	partition
        0	1 | 2 3 4 
        1	2 | 1 3 4 
        2	3 | 1 2 4 
        3	4 | 1 2 3 
        4	1 2 | 3 4 
"""
function show_bipartitions(n::Int64; start::Int64 = 0, stop::Int64=-1)
    idx_width = length(string(num_bipartitions(n)))
    idx_fmt = FormatSpec(string(">", idx_width, "s"))
    node_width = length(string(n))
    node_fmt = FormatSpec(string(">", node_width, "s"))
    idx = 0
    node = Vector(1:n)
    println("idx\tpartition")
    for i in range(1, n÷2)
        
        comb = combinations(Vector(1:n), i)
        if i == n - i
            comb = Iterators.take(comb,  binomial(n,i)÷ 2)
        end
            
        for c in comb
            if start <= idx && (stop == -1 || idx <= stop)
                printfmt(idx_fmt, idx)
                print("\t")
                for e in c
                    printfmt(node_fmt, e)
                    print(" ")
                end
                print("| ")
                for e in sort(collect(setdiff(node ,Set(c))))
                    printfmt(node_fmt, e)
                    print(" ")
                end
                println()
            end
            idx += 1
        end
    end
end;

"""
        show_bipartition(n::Int64, idx::Int64)

    Function to show the bipartition with a given index of the number of taxa.

    # Arguments
        n: the number of taxa
        idx: the bipartition to show
    # Output
        the idx-th bipartiion of trees with n taxa
        Example:
        show_bipartition(4, 3):
            idx	partition
            3	4 | 1 2 3 
"""
function show_bipartition(n::Int64, idx::Int64)
    show_bipartitions(n,start = idx,stop = idx)
end;

"""
    Function to get existing taxa on an edge
    Input:
        encoded_taxa: a Vector of booleans where true for taxa that are descendent of the edge, false for other taxa
    output:
        a Vector of taxa that exist on the edge
"""
function get_taxa(encoded_taxa)
    nodes = []
    for i in 1:length(encoded_taxa)
        if encoded_taxa[i]
            append!(nodes, i)
        end
    end
    return nodes
end

"""
    Function to get the bipartition format of a phylogenetic tree
    Input:
        tree: a HybridNetwork object that reperents the tree
        n: the number of taxa of the tree
    output:
        a Vector of pairs shows the tree in bipartiton format
    Example:
        tree:
            HybridNetwork, Rooted Network
            6 edges
            7 nodes: 4 tips, 0 hybrid nodes, 3 internal tree nodes.
            tip labels: 4, 1, 2, 3
            (4:4.249,(1:2.457,(2:2.064,3:2.064):0.393):1.793);

        get_bipartition(tree, 4):
            6-element Vector{Any}:
                3 => 4.249
                0 => 2.457
                1 => 2.064
                2 => 2.064
                6 => 0.393
                3 => 1.793
"""
function get_bipartition(tree::HybridNetwork, n::Int64)
    taxa = sort(tipLabels(tree))
    if length(taxa) != n
        print("check the parameter n")
        return
    end
    node = Vector(1:n)
    result = []
    idx = 0
    for i in tree.edge
        node_idx = hardwiredCluster(i,taxa)
        branch_node = get_taxa(node_idx)
        
        # if branch seperates more than half nodes, we use the small part to get bipartition idx
        if (length(branch_node) > n ÷ 2)
             branch_node = sort(collect(setdiff(node,branch_node)))
        end
        
        # generate all possible combination with the same number of nodes
        comb = collect(combinations(Vector(1:n), length(branch_node)))
        for c in 1:length(comb)
            if comb[c] == branch_node
                # if the combination is the later, we need to find its first half
                if length(branch_node) > 1 && c > length(comb) ÷ 2
                    c = length(comb) - c + 1 
                end
                idx = c - 1
                break
            end
        end
        for j in 1:(length(branch_node) - 1)
            idx = idx + binomial(n,j)
        end
        push!(result,Pair(idx, i.length))
        idx = 0
    end
    return result
end


"""
    Function to get the bipartition format of phylogenetic trees in Matrix
    Note that all taxa will be replaced by numbers. Use the function num_to_name to 
    get a Dictionary containing the name of the taxon corresponding to the number
    Input:
        trees: a Vector of HybridNetwork objects that contains trees
        n: the number of taxa of the trees
    output:
        a Matrix shows the tree in bipartiton format
    Example:
        tree:
            HybridNetwork, Rooted Network
            6 edges
            7 nodes: 4 tips, 0 hybrid nodes, 3 internal tree nodes.
            tip labels: 4, 1, 2, 3
            (4:4.249,(1:2.457,(2:2.064,3:2.064):0.393):1.793);

        print_bipartition([tree], 4)
            1*7 Matrix{Float64}:
            2.457  2.064  2.064  6.042  0.0  0.0  0.393
"""
function print_bipartition(trees::Vector{HybridNetwork}, n::Int64)
    
    # build the table
    N = num_bipartitions(n)
    data = zeros(length(trees), N)
    treeNum = 1;
    # get existing index
    
    for tree in trees
        bipart = get_bipartition(tree, n)
        for j in 1:length(bipart)
            data[treeNum,(bipart[j][1] + 1)] += bipart[j][2]
        end
        treeNum+=1
    end
    return data
end

"""
    Function to obtain a table that contains taxon numbers and their names 
    Input:
        trees: a HybridNetwork objects that contains trees
    output:
        a Dictionary shows taxon numbers and their names
    Example:
        tree:
            HybridNetwork, Rooted Network
            6 edges
            7 nodes: 4 tips, 0 hybrid nodes, 3 internal tree nodes.
            tip labels: O, HYB, P1, P2
            (O:3.866,(HYB:1.593,(P1:1.208,P2:1.208):0.386):2.273);

        print_bipartition([tree], 4)
            Dict{Int64, Any} with 4 entries:
            4 => "P2"
            2 => "O"
            3 => "P1"
            1 => "HYB"
"""
function num_to_name(tree::HybridNetwork)
    taxa = sort(tipLabels(tree))
    dict = Dict{Int64, Any}()
    for i in 1:length(taxa)
        dict[i] = taxa[i]
    end
    return dict
end