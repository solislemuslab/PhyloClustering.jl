using Clustering, GaussianMixtures, ParallelKMeans, Random

"""
    kmeans_label(tree::AbstractMatrix{<:Real}, n::Int64; init::String="k-means++", rng::AbstractRNG=Random.GLOBAL_RNG)

Function to get predicted labels from Yinyang K-means clustering.

# Arguments
 - `tree`: a B * N tree Matrix (each column of tree Matrix is a B-dimensional tree in bipartiton format).
 - `n`: the number of clusters.
 - `init`(defaults to `k-means++`): is one of the algorithms used for initialization.
    Alternatively, one can use `rand` to choose random points for init.
 - `rng`: RNG object.
# Output
 A `Vector` object with length of N containing predicted labels for each tree (the cluster it belongs to).
"""
function kmeans_label(tree::AbstractMatrix{<:Real}, n::Int64; init::String="k-means++", rng::AbstractRNG=Random.GLOBAL_RNG)  
    result = ParallelKMeans.kmeans(Yinyang(), tree, n; k_init =init, rng=rng)
    return result.assignments
end

"""
    gmm_label(tree::AbstractMatrix{<:Real}, n::Int64; method::Symbol=:kmeans, kind::Symbol=:diag)

Function to get predicted labels from Gaussian mixture model.

# Arguments
 - `tree`: a B * N tree Matrix (each column of tree Matrix is a B-dimensional tree in bipartiton format).
 - `n`: the number of clusters.
 - `method`(defaults to `:kmeans`): intialization method to find n starting centers:
        * `:kmeans`: use K-means clustering from [Clustering.jl](https://github.com/JuliaStats/Clustering.jl) to initialize with `n` centers.
        * `:split`: initialize a single Gaussian with `tree` and subsequently splitting the Gaussians followed by retraining using 
        the EM algorithm until `n` Gaussians are obtained. 
 - `kind`(defaults to `:diag`): covariance type, `:diag` or `:full`.
# Output
 A `Tuple{Vector{Int64}, Vector{Int64}}` where the first `Vector` contains predicted labels for each tree based on the posterior probability and
    the second `Vector` contain predicted labels for each tree based on the Log Likelihood.
"""
function gmm_label(tree::AbstractMatrix{<:Real}, n::Int64; method::Symbol=:kmeans, kind::Symbol=:diag)
    data= Matrix(tree');
    gmm= GMM(n,data,method=method, kind=kind);
    prob_pos=gmmposterior(gmm,data)[1]
    llpg = gmmposterior(gmm,data)[2]
    llpg_result = map(argmax, eachrow(llpg))
    prob_result = map(argmax, eachrow(prob_pos))
    return prob_result, llpg_result
end;

"""
    hc_label(matrix::AbstractMatrix{<:Real}, n::Int64; linkage::Symbol=:ward)

Function to get predicted labels from hierarchical clustering.

# Arguments
 - `matrix`: a N * N pairwise distance Matrix.
 - `n`: the number of clusters.
 - `linkage`(defaults to `:ward`): cluster linkage function to use. It affects what clusters are merged on each iteration:
   * `:single`: use the minimum distance between any of the cluster members
   * `:average`: use the mean distance between any of the cluster members
   * `:complete`: use the maximum distance between any of the members
   * `:ward`: the distance is the increase of the average squared distance of a point to its cluster centroid after merging the two clusters
   * `:ward_presquared`: same as `:ward`, but assumes that the distances in `d` are already squared.
# Output
 A `Vector` object with length of N containing predicted labels for each tree (the cluster it belongs to). 
"""
function hc_label(matrix::AbstractMatrix{<:Real}, n::Int64; linkage::Symbol=:ward)
    H = hclust(matrix,linkage = linkage)
    return pred = cutree(H, k=n)
end

"""
    dbscan_label(tree::AbstractMatrix{<:Real}, radius::Real; min_neighbors::Int64 = 1, min_cluster_size::Int64 = 1)
    
Function to get predicted labels from density-based spatial clustering of applications with noise.

# Arguments
 - `tree`: a B * N tree Matrix (each column of tree Matrix is a B-dimensional tree in bipartiton format) and B < N.
 - `radius`: neighborhood radius; points within this distance are considered neighbors.
 - `min_neighbors`: minimal number of neighbors required to assign a point to a cluster.
 - `min_cluster_size`: minimal number of points in a cluster.
# Output
 A `Vector` object with length of N containing predicted labels for each tree (the cluster it belongs to). 
    0 means the tree is noise.
"""
function dbscan_label(tree::AbstractMatrix{<:Real}, radius::Real; min_neighbors::Int64 = 1, min_cluster_size::Int64 = 1)   
    result = dbscan(tree, radius,min_neighbors = min_neighbors, min_cluster_size = min_cluster_size)   
    # get only points in clusters
    return result.assignments
end
