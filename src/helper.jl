using MultivariateStats, StatsBase, Distances

"""
    standardize_tree(tree::AbstractMatrix{<:Real})

Standardize tree Matrix that returned by [`split_weight`](@ref). 
It is recommended to standardize the data before inputting it into the model.

# Arguments
 - `tree`: a N * B Matrix containing trees (each row is a B-dimensional tree in bipartiton format).
# Output
 A standardized B * N tree `Matrix` with a mean of about 0 and a standard deviation of about 1.
    This tree `Matrix` can be the input of [model](@ref basics).
"""
function standardize_tree(tree::AbstractMatrix{<:Real})
    data = collect(tree')
    dt = fit(ZScoreTransform, data, dims=2)
    data = StatsBase.transform(dt, data)
    replace!(data, NaN => 0)
    return data
end

"""
    distance(tree::AbstractMatrix{<:Real})
    
Get the distance Matrix of a tree Matrix returned by [`split_weight`](@ref).

# Arguments
 - `tree`: a B * N tree Matrix (each column of tree Matrix is a B-dimensional tree in bipartiton format).
# Output
 A pairwise distance `Matrix` that can be the input of [`hc_label`](@ref).
"""
function distance(tree::AbstractMatrix{<:Real})
    matrix = pairwise(Euclidean(), tree, dims=2)
    return matrix
end