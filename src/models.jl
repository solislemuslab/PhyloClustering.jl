using Clustering, GaussianMixtures, ParallelKMeans

"""
    Function to get predicted labels from Yinyang K-means clustering
    Input:
        tree: a B * N tree Matrix (each column of tree Matrix is a B-dimensional tree in bipartiton format)
        n: the number of clusters
        init: initialization method. K-means++ or rand
    output:
        a Vector object with length of N containing predicted labels for each tree
"""
function kmeans_label(tree::AbstractMatrix{<:Real}, n::Int64; init::String="k-means++")  
    result = ParallelKMeans.kmeans(Yinyang(), tree, n; k_init =init)
    return result.assignments
end

"""
    Function to get predicted labels from Gaussian mixture model.
    Input:
        tree: a B * N tree Matrix (each column of tree Matrix is a B-dimensional tree in bipartiton format)
        n: the number of clusters
        method: intialization method to find n starting centers
        kind: covariance type
    output:
        a vector tuple where the first vector contains predicted labels for each tree based on the posterior probability and
        the second vector contain predicted labels for each tree based on the Log Likelihood 
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
    Function to get predicted labels from hierarchical clustering.
    Input:
        matrix: a N * N pairwise distance Matrix
        n: the number of clusters
        linkage: cluster method to use. It affects what clusters are merged on each iteration.
    output:
        a Vector object with length of N containing predicted labels for each tree (the cluster it belongs to) 
"""
function hc_label(matrix::AbstractMatrix{<:Real}, n::Int64; linkage::Symbol=:ward)
    H = hclust(matrix,linkage = linkage)
    return pred = cutree(H, k=n)
end

"""
    Function to get predicted labels from hierarchical clustering.
    Input:
        tree: a B * N tree Matrix (each column of tree Matrix is a B-dimensional tree in bipartiton format)
        radius: neighborhood radius; points within this distance are considered neighbors
        min_neighbors: minimal number of neighbors required to assign a point to a cluster
        min_cluster_size: minimal number of points in a cluster
    output:
        a Vector object with length of N containing predicted labels for each tree (the cluster it belongs to) 
"""
function dbscan_label(tree::AbstractMatrix{<:Real}, radius::Real; min_neighbors::Int64 = 1, min_cluster_size::Int64 = 1)   
    result = dbscan(tree, radius,min_neighbors = min_neighbors, min_cluster_size = min_cluster_size)   
    # get only points in clusters
    result = getproperty.(result, :core_indices)
    idx = fill(0,length(tree[1,:]))
    for i in range(1, length(idx))
        if i in result[1]
            idx[i] = 1
        elseif i in result[2]
            idx[i] = 2
        end
    end
    return idx
end
