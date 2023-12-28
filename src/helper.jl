using MultivariateStats, StatsBase, CairoMakie, Distances

"""
    standardize_tree(tree::AbstractMatrix{<:Real})

Function to standardize tree Matrix.

# Arguments
 - `tree`: a N * B Matrix containing trees (each row is a B-dimensional tree in bipartiton format).
# Output
 A standardized B * N tree `Matrix` with a mean of about 0 and a standard deviation of about 1.
    This tree `Matrix` can be the input of [model](@ref basics).
"""
function standardize_tree(tree::AbstractMatrix{<:Real})    
    data = collect(tree');
    dt = fit(ZScoreTransform, data, dims=2)
    data = StatsBase.transform(dt, data)
    replace!(data, NaN=>0)
    return data
end

"""
    plot_clusters(tree::AbstractMatrix{<:Real}, label::Vector{Int64})

Function to visualize the result of models.

# Arguments
 - `tree`: a B * N tree Matrix (each column of tree Matrix is a B-dimensional tree in bipartiton format).
 - `label`: an N-length Vector containing predicted labels for each tree. People can use the output of the models.
# Output
 A scatter plot showing tree clusters.
"""
function plot_clusters(tree::AbstractMatrix{<:Real}, label::Vector{Int64})
    PCA_model = fit(PCA, tree, maxoutdim = 2);
    PCA_data = predict(PCA_model,tree)
    scatter(PCA_data[1,:], PCA_data[2,:], markersize = 5, color = label)
end

"""
    distance(tree::AbstractMatrix{<:Real})
    
Function to get the distance Matrix of a tree Matrix.

# Arguments
 - `tree`: a B * N tree Matrix (each column of tree Matrix is a B-dimensional tree in bipartiton format).
# Output
 A pairwise distance `Matrix` that can be the input of [`hc_label`](@ref).
"""
function distance(tree::AbstractMatrix{<:Real})
    matrix = pairwise(Euclidean(), tree, dims=2)
    return matrix
end