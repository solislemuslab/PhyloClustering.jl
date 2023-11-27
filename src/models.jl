using MultivariateStats, StatsBase, GaussianMixtures, ParallelKMeans, MLBase, Hungarian, LinearAlgebra

"""
    Function to get predicted labels from Yinyang K-means clustering
    Input:
        tree: a B * N tree matrix (each column of tree matrix is a B-dimensional tree in bipartiton format)
        n: the number of clusters
        init: initialization method. K-means++ or rand
    output:
        a Vector object with length of N containing the predicted label of each tree (the cluster it belongs to)
"""
function kmeans_label(tree, n; init=:"k-means++")  
    result = ParallelKMeans.kmeans(Yinyang(), tree, n; k_init =init)
    return result.assignments
end

"""
    Function to get predicted labels from Gaussian mixture model.
    Input:
        tree: a B * N tree matrix (each column of tree matrix is a B-dimensional tree in bipartiton format)
        n: the number of clusters
        method: intialization method to find n starting centers
        kind: covariance type
    output:
        a Vector object with length of N containing the predicted label of each tree (the cluster it belongs to)
"""
function GMM_label(tree, n; method=:kmeans, kind=:diag)    
    data= tree'
    gmm=GMM(n,Array(data),method=method, kind=kind);
    prob_pos=gmmposterior(gmm,Array(data))[1]
    llpg = gmmposterior(gmm,Array(data))[2]
    llpg_result = map(argmax, eachrow(llpg))
    prob_result = map(argmax, eachrow(prob_pos))
    return prob_result, llpg_result
end

function hc_label(matrix, n; linkage=:ward)
    H = hclust(matrix,linkage = linkage)
    return pred = cutree(H, k=n)
end

function dbscan_label(tree, radius; min_neighbors = 1, min_cluster_size = 1)   
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
