using MultivariateStats, StatsBase, GaussianMixtures, ParallelKMeans, MLBase, Hungarian, LinearAlgebra

"""
    Function to get predicted labels from K-means clustering
    Input:
        tree: DataFrame object with the concordance factor table (first 4 columns should be taxon names, last 3 columns should be cf values)
        n: the number of clusters
        seed: 
    output:
        top m optimal phylogenetic networks in a Dict object
"""
function kmeans_label(tree, n; seed =:"k-means++")  
    result = ParallelKMeans.kmeans(Yinyang(),tree, n; k_init ="k-means++")
    return result.assignments
end

function GMM_label(tree, n; method=:kmeans, kind=:diag)    
    data= tree'
    gmm=GMM(n,Array(data),method=method, kind=kind);
    prob_pos=gmmposterior(gmm,Array(data))[1]
    llpg = gmmposterior(gmm,Array(data))[2]
    llpg_result = map(argmax, eachrow(llpg))
    prob_result = map(argmax, eachrow(prob_pos))
    return Pair(llpg_result, prob_result)
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
