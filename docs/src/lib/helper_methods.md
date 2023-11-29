```@meta
CurrentModule = PhyloClustering
```

# Convert trees in Newick format to bipartition format
```@docs
num_bipartitions(n::Int64)
```

```@docs
show_bipartitions(n::Int64; start::Int64 = 0, stop::Int64=-1)
```

```@docs
show_bipartition(n::Int64, idx::Int64)
```

```@docs
plot_clusters(tree::AbstractMatrix{<:Real}, label::Vector{Int64})
```

```@docs
get_bipartition(tree::HybridNetwork, n::Int64)
```

```@docs
print_bipartition(trees::Vector{HybridNetwork}, n::Int64)
```

```@docs
standardize_tree(tree::AbstractMatrix{<:Real}) 
```

```@docs
num_to_name(tree::HybridNetwork)
```

```@docs
distance(tree::AbstractMatrix{<:Real})
```