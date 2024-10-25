var documenterSearchIndex = {"docs":
[{"location":"model/dbscan/#DBSCAN","page":"DBSCAN","title":"DBSCAN","text":"","category":"section"},{"location":"model/dbscan/","page":"DBSCAN","title":"DBSCAN","text":"Density-based Spatial Clustering of Applications with Noise (DBSCAN) is a data clustering algorithm that finds clusters through density-based expansion of seed points. The algorithm was proposed in:","category":"page"},{"location":"model/dbscan/","page":"DBSCAN","title":"DBSCAN","text":"dbscan_label","category":"page"},{"location":"model/dbscan/#PhyloClustering.dbscan_label","page":"DBSCAN","title":"PhyloClustering.dbscan_label","text":"dbscan_label(tree::AbstractMatrix{<:Real}, radius::Real; min_neighbors::Int64 = 1, min_cluster_size::Int64 = 1)\n\nGet predicted labels from density-based spatial clustering of applications with noise for a group of phylogenetic trees.\n\nArguments\n\ntree: a B * N tree Matrix (each column of tree Matrix is a B-dimensional tree in bipartiton format) and B < N.\nradius: neighborhood radius; points within this distance are considered neighbors.\nmin_neighbors: minimal number of neighbors required to assign a point to a cluster.\nmin_cluster_size: minimal number of points in a cluster.\n\nOutput\n\nA Vector object with length of N containing predicted labels for each tree (the cluster it belongs to).      0 means the tree is noise.\n\n\n\n\n\n","category":"function"},{"location":"model/dbscan/","page":"DBSCAN","title":"DBSCAN","text":"Reference","category":"page"},{"location":"model/dbscan/","page":"DBSCAN","title":"DBSCAN","text":"The implementation of DBSCAN is provided by Clustering.jl.","category":"page"},{"location":"model/dbscan/","page":"DBSCAN","title":"DBSCAN","text":"Martin Ester, Hans-Peter Kriegel, Jörg Sander, and Xiaowei Xu. 1996.  A density-based algorithm for discovering clusters in large spatial databases with noise.  In Proceedings of the Second International Conference on Knowledge Discovery and Data Mining (KDD'96). AAAI Press, 226–231.","category":"page"},{"location":"model/basic/#basics","page":"Basics","title":"Basics","text":"","category":"section"},{"location":"model/basic/","page":"Basics","title":"Basics","text":"The package implements a variety of clustering algorithms:","category":"page"},{"location":"model/basic/","page":"Basics","title":"Basics","text":"Pages = [\"kmeans.md\", \"hclust.md\", \"dbscan.md\", \"gmm.md\"]","category":"page"},{"location":"model/basic/","page":"Basics","title":"Basics","text":"K-means, DBSCAN, and GMM need the same input tree matrix, making it easy to use different models. Check example for the usage of these three models.","category":"page"},{"location":"lib/bipartition/#Split-Weight-Embedding-of-Phylogenetic-trees","page":"Split-Weight Embedding","title":"Split-Weight Embedding of Phylogenetic trees","text":"","category":"section"},{"location":"lib/bipartition/","page":"Split-Weight Embedding","title":"Split-Weight Embedding","text":"Our package provides functions to convert phylogenetic trees from Newick format to Bipartition format, which are then embedded in Euclidean space via Split-Weight Embedding. ","category":"page"},{"location":"lib/bipartition/","page":"Split-Weight Embedding","title":"Split-Weight Embedding","text":"Check our paper for more details.","category":"page"},{"location":"lib/bipartition/","page":"Split-Weight Embedding","title":"Split-Weight Embedding","text":"num_bipartitions\nshow_bipartitions\nshow_bipartition\nget_bipartition\nsplit_weight\nnum_to_name","category":"page"},{"location":"lib/bipartition/#PhyloClustering.num_bipartitions","page":"Split-Weight Embedding","title":"PhyloClustering.num_bipartitions","text":"num_bipartitions(n::Int64)\n\nGet the number of bipartitions for a given number of taxa.\n\nArguments\n\nn: the number of taxa.\n\nOutput:\n\nThe number of bipartition corresponding to the number of taxa.\n\n\n\n\n\n","category":"function"},{"location":"lib/bipartition/#PhyloClustering.show_bipartitions","page":"Split-Weight Embedding","title":"PhyloClustering.show_bipartitions","text":"show_bipartitions(n::Int64; start::Int64 = 0, stop::Int64=-1)\n\nShow the bipartitions between certain indices for a given number of taxa.\n\nArguments\n\nn: the number of taxa\nstart(defaults to 0): the starting index.\nstop(defaults to -1): the ending index.\n\nOutput\n\nThe bipartiions of trees with n taxa between start index and stop index.\n\nExample\n\njulia> show_bipartitions(4,stop = 4)\nidx\tpartition\n0\t1 | 2 3 4 \n1\t2 | 1 3 4 \n2\t3 | 1 2 4 \n3\t4 | 1 2 3 \n4\t1 2 | 3 4 \n\n\n\n\n\n","category":"function"},{"location":"lib/bipartition/#PhyloClustering.show_bipartition","page":"Split-Weight Embedding","title":"PhyloClustering.show_bipartition","text":"show_bipartition(n::Int64, idx::Int64)\n\nShow the bipartition with a given index for a given number of taxa.\n\nArguments\n\nn: the number of taxa.\nidx: the bipartition to show.\n\nOutput\n\nThe idx-th bipartiion of trees with n taxa.\n\nExample:\n\njulia> show_bipartition(4, 3)\nidx\tpartition\n3\t4 | 1 2 3 \n\n\n\n\n\n","category":"function"},{"location":"lib/bipartition/#PhyloClustering.get_bipartition","page":"Split-Weight Embedding","title":"PhyloClustering.get_bipartition","text":"get_bipartition(tree::HybridNetwork, n::Int64)\n\nGet the bipartition format of a phylogenetic tree.\n\nArguments\n\ntree: a HybridNetwork object that reperents the tree.\nn: the number of taxa of the tree.\n\nOutput\n\nA Vector of pairs shows the tree in bipartiton format.\n\nExample\n\njulia> get_bipartition(readTopology(\"(4:4.249,(1:2.457,(2:2.064,3:2.064):0.393):1.793);\"), 4)\n6-element Vector{Any}:\n 3 => 4.249\n 0 => 2.457\n 1 => 2.064\n 2 => 2.064\n 6 => 0.393\n 3 => 1.793\n\n\n\n\n\n","category":"function"},{"location":"lib/bipartition/#PhyloClustering.split_weight","page":"Split-Weight Embedding","title":"PhyloClustering.split_weight","text":"split_weight(trees::Vector{HybridNetwork}, n::Int64)\n\nSplit-weight embedding of phylogenetic trees in a Vector{HybridNetwork}. Note that all taxa will be replaced by numbers. Use the function num_to_name to  get a Dictionary containing the name of the taxon corresponding to the number.\n\nArguments\n\ntrees: a Vector of HybridNetwork objects that contains trees\nn: the number of taxa of the trees\n\nOutput\n\nA Matrix{Float64} shows the trees in bipartiton format\n\nExample\n\njulia> split_weight([readTopology(\"(4:4.249,(1:2.457,(2:2.064,3:2.064):0.393):1.793);\"), readTopology(\"(4:4.308,(2:1.634,(1:1.588,3:1.588):0.046):2.674);\")], 4)\n2×7 Matrix{Float64}:\n 2.457  2.064  2.064  6.042  0.0  0.0    0.393\n 1.588  1.634  1.588  6.982  0.0  0.046  0.0\n\n\n\n\n\n","category":"function"},{"location":"lib/bipartition/#PhyloClustering.num_to_name","page":"Split-Weight Embedding","title":"PhyloClustering.num_to_name","text":"num_to_name(tree::HybridNetwork)\n\nObtain a table that contains taxon numbers and their names for a given tree.\n\nArguments\n\ntree: a HybridNetwork objects that contains trees\n\nOutput\n\nA Dict{Int64, Any} shows taxon numbers and their names\n\nExample\n\njulia> num_to_name(readTopology(\"(O:3.866,(HYB:1.593,(P1:1.208,P2:1.208):0.386):2.273);\"))\nDict{Int64, Any} with 4 entries:\n  4 => \"P2\"\n  2 => \"O\"\n  3 => \"P1\"\n  1 => \"HYB\"\n\n\n\n\n\n","category":"function"},{"location":"model/gmm/#Gaussian-Mixture-Model-(GMM)","page":"Gaussian Mixture Model (GMM)","title":"Gaussian Mixture Model (GMM)","text":"","category":"section"},{"location":"model/gmm/","page":"Gaussian Mixture Model (GMM)","title":"Gaussian Mixture Model (GMM)","text":"Gaussian Mixture Model (GMM) a model-based probabilistic method that assumes data points are generated from a mixture of several Gaussian distributions with unknown parameters. It uses the Expectation Maximization (EM) algorithm to update the parameters iteratively in order to optimize the log-likelihood of the data until convergence.","category":"page"},{"location":"model/gmm/","page":"Gaussian Mixture Model (GMM)","title":"Gaussian Mixture Model (GMM)","text":"gmm_label","category":"page"},{"location":"model/gmm/#PhyloClustering.gmm_label","page":"Gaussian Mixture Model (GMM)","title":"PhyloClustering.gmm_label","text":"gmm_label(tree::AbstractMatrix{<:Real}, n::Int64; method::Symbol=:kmeans, kind::Symbol=:diag)\n\nGet predicted labels from Gaussian mixture model for a group of phylogenetic trees.\n\nArguments\n\ntree: a B * N tree Matrix (each column of tree Matrix is a B-dimensional tree in bipartiton format).\nn: the number of clusters.\nmethod(defaults to :kmeans): intialization method to find n starting centers:      * :kmeans: use K-means clustering from Clustering.jl to initialize with n centers.      * :split: initialize a single Gaussian with tree and subsequently splitting the Gaussians followed by retraining using       the EM algorithm until n Gaussians are obtained. \nkind(defaults to :diag): covariance type, :diag or :full.\n\nOutput\n\nA Tuple{Vector{Int64}, Vector{Int64}} where the first Vector contains predicted labels for each tree based on the posterior probability and     the second Vector contain predicted labels for each tree based on the Log Likelihood.\n\n\n\n\n\n","category":"function"},{"location":"model/gmm/","page":"Gaussian Mixture Model (GMM)","title":"Gaussian Mixture Model (GMM)","text":"Reference","category":"page"},{"location":"model/gmm/","page":"Gaussian Mixture Model (GMM)","title":"Gaussian Mixture Model (GMM)","text":"The implementation of GMM is provided by GaussianMixtures.jl.","category":"page"},{"location":"model/hclust/#Hierarchical-Clustering","page":"Hierarchical Clustering","title":"Hierarchical Clustering","text":"","category":"section"},{"location":"model/hclust/","page":"Hierarchical Clustering","title":"Hierarchical Clustering","text":"Hierarchical clustering is a hierarchical method that creates a clustering tree called a dendrogram. Each leaf on the tree is a data point, and branches represent the joining of clusters. We can 'cut' the tree at different heights to get a different number of clusters. This method does not require the number of clusters to be specified in advance and can be either agglomerative (bottom-up) or divisive (top-down).","category":"page"},{"location":"model/hclust/","page":"Hierarchical Clustering","title":"Hierarchical Clustering","text":"hc_label","category":"page"},{"location":"model/hclust/#PhyloClustering.hc_label","page":"Hierarchical Clustering","title":"PhyloClustering.hc_label","text":"hc_label(matrix::AbstractMatrix{<:Real}, n::Int64; linkage::Symbol=:ward)\n\nGet predicted labels from hierarchical clustering for a group of phylogenetic trees.\n\nArguments\n\nmatrix: a N * N pairwise distance Matrix.\nn: the number of clusters.\nlinkage(defaults to :ward): cluster linkage function to use. It affects what clusters are merged on each iteration:\n:single: use the minimum distance between any of the cluster members\n:average: use the mean distance between any of the cluster members\n:complete: use the maximum distance between any of the members\n:ward: the distance is the increase of the average squared distance of a point to its cluster centroid after merging the two clusters\n:ward_presquared: same as :ward, but assumes that the distances in d are already squared.\n\nOutput\n\nA Vector object with length of N containing predicted labels for each tree (the cluster it belongs to). \n\n\n\n\n\n","category":"function"},{"location":"model/hclust/#Example","page":"Hierarchical Clustering","title":"Example","text":"","category":"section"},{"location":"model/hclust/","page":"Hierarchical Clustering","title":"Hierarchical Clustering","text":"using PhyloClustering, PhyloNetworks\n\n# read trees with 4-taxa in Newick format using PhyloNetworks\ntrees = readMultiTopology(\"../data/data.trees\");\n\n# convert trees to Bipartition foramt and embed them via split-weight embedding\ntrees = split_weight(trees, 4);","category":"page"},{"location":"model/hclust/","page":"Hierarchical Clustering","title":"Hierarchical Clustering","text":"Standardize the data, calculate the distance matrix, and input them into Yinyang K-means clustering.","category":"page"},{"location":"model/hclust/","page":"Hierarchical Clustering","title":"Hierarchical Clustering","text":"tree = standardize_tree(trees);\nmatrix = distance(tree);\nlabel = hc_label(matrix, 2)","category":"page"},{"location":"model/hclust/","page":"Hierarchical Clustering","title":"Hierarchical Clustering","text":"Reference","category":"page"},{"location":"model/hclust/","page":"Hierarchical Clustering","title":"Hierarchical Clustering","text":"The implementation of hierarchical clustering is provided by Clustering.jl.","category":"page"},{"location":"model/hclust/","page":"Hierarchical Clustering","title":"Hierarchical Clustering","text":"Joe H. Ward Jr. (1963) Hierarchical Grouping to Optimize an Objective Function,  Journal of the American Statistical Association, 58:301, 236-244, DOI: 10.1080/01621459.1963.10500845","category":"page"},{"location":"model/kmeans/#K-means","page":"K-means","title":"K-means","text":"","category":"section"},{"location":"model/kmeans/","page":"K-means","title":"K-means","text":"K-means is a classical method for clustering or vector quantization. It produces a fixed number of clusters, each associated with a center (also known as a prototype), and each data point is assigned to a cluster with the nearest center.","category":"page"},{"location":"model/kmeans/","page":"K-means","title":"K-means","text":"The K-means clustering used in our package is Yinyang K-means which has less runtime and memory usage on large datasets. Traditional K-means calculates the distance from all data points to centroids for each iteration. Yinyang K-means uses triangular inequalities to construct and maintain upper and lower bounds on the distances of data points from the centroids, with global and local filtering to minimize unneeded calculations.","category":"page"},{"location":"model/kmeans/","page":"K-means","title":"K-means","text":"kmeans_label","category":"page"},{"location":"model/kmeans/#PhyloClustering.kmeans_label","page":"K-means","title":"PhyloClustering.kmeans_label","text":"kmeans_label(tree::AbstractMatrix{<:Real}, n::Int64; init::String=\"k-means++\", rng::AbstractRNG=Random.GLOBAL_RNG)\n\nGet predicted labels from Yinyang K-means clustering for a group of phylogenetic trees.\n\nArguments\n\ntree: a B * N tree Matrix (each column of tree Matrix is a B-dimensional tree in bipartiton format).\nn: the number of clusters.\ninit(defaults to k-means++): is one of the algorithms used for initialization.  Alternatively, one can use rand to choose random points for init.\nrng: RNG object.\n\nOutput\n\nA Vector object with length of N containing predicted labels for each tree (the cluster it belongs to).\n\n\n\n\n\n","category":"function"},{"location":"model/kmeans/#usage","page":"K-means","title":"Example","text":"","category":"section"},{"location":"model/kmeans/","page":"K-means","title":"K-means","text":"using PhyloClustering, PhyloNetworks, StableRNGs\n# use RNG for stable result\nrng = StableRNG(1)\n\n# read trees with 4-taxa in Newick format using PhyloNetworks\ntrees = readMultiTopology(\"../data/data.trees\");\n\n# convert trees to Bipartition foramt and embed them via split-weight embedding\ntrees = split_weight(trees, 4)","category":"page"},{"location":"model/kmeans/","page":"K-means","title":"K-means","text":"Standardize the data and input them into Yinyang K-means clustering.","category":"page"},{"location":"model/kmeans/","page":"K-means","title":"K-means","text":"tree = standardize_tree(trees);\nlabel = kmeans_label(tree, 2, rng=rng)","category":"page"},{"location":"model/kmeans/","page":"K-means","title":"K-means","text":"Reference","category":"page"},{"location":"model/kmeans/","page":"K-means","title":"K-means","text":"The implementation of Yinyang K-means is provided by ParallelKMeans.jl.","category":"page"},{"location":"model/kmeans/","page":"K-means","title":"K-means","text":"Yufei Ding et al. 2015. Yinyang K-Means: A Drop-In Replacement of the Classic K-Means with Consistent Speedup.  Proceedings of the 32nd International Conference on Machine Learning, ICML 2015, Lille, France, 6-11 July 2015","category":"page"},{"location":"#PhyloClustering.jl","page":"Home","title":"PhyloClustering.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PhyloClustering.jl is a Julia package for performing unsupervised learning on phylogenetic trees. The algorithms currently included are K-means, Hierarchical Clustering, Gaussian Mixture Model (GMM), and Density-based Spatial Clustering of Applications with Noise (DBSCAN).","category":"page"},{"location":"#Citation","page":"Home","title":"Citation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you use PhyloClustering.jl in your work, we kindly ask that you cite the following paper: ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Kong, Y., Tiley, G. P., Solís-Lemus, C. (2023). Unsupervised learning of phylogenetic trees via split-weight embedding. arXiv:2312.16074.","category":"page"},{"location":"lib/helper_methods/#Pre-process-data-before-input-models","page":"Helper Functions","title":"Pre-process data before input models","text":"","category":"section"},{"location":"lib/helper_methods/","page":"Helper Functions","title":"Helper Functions","text":"standardize_tree\ndistance","category":"page"},{"location":"lib/helper_methods/#PhyloClustering.standardize_tree","page":"Helper Functions","title":"PhyloClustering.standardize_tree","text":"standardize_tree(tree::AbstractMatrix{<:Real})\n\nStandardize tree Matrix that returned by split_weight.  It is recommended to standardize the data before inputting it into the model.\n\nArguments\n\ntree: a N * B Matrix containing trees (each row is a B-dimensional tree in bipartiton format).\n\nOutput\n\nA standardized B * N tree Matrix with a mean of about 0 and a standard deviation of about 1.     This tree Matrix can be the input of model.\n\n\n\n\n\n","category":"function"},{"location":"lib/helper_methods/#PhyloClustering.distance","page":"Helper Functions","title":"PhyloClustering.distance","text":"distance(tree::AbstractMatrix{<:Real})\n\nGet the distance Matrix of a tree Matrix returned by split_weight.\n\nArguments\n\ntree: a B * N tree Matrix (each column of tree Matrix is a B-dimensional tree in bipartiton format).\n\nOutput\n\nA pairwise distance Matrix that can be the input of hc_label.\n\n\n\n\n\n","category":"function"},{"location":"man/installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"man/installation/#Installation-of-Julia","page":"Installation","title":"Installation of Julia","text":"","category":"section"},{"location":"man/installation/","page":"Installation","title":"Installation","text":"Julia is a high-level and interactive programming language (like R or Matlab), but it is also high-performance (like C). To install Julia, follow instructions here. For a quick & basic tutorial on Julia, see learn x in y minutes.","category":"page"},{"location":"man/installation/","page":"Installation","title":"Installation","text":"Editors:","category":"page"},{"location":"man/installation/","page":"Installation","title":"Installation","text":"Visual Studio Code provides an editor and an integrated development environment (IDE) for Julia: highly recommended!\nYou can also run Julia within a Jupyter notebook (formerly IPython notebook).","category":"page"},{"location":"man/installation/","page":"Installation","title":"Installation","text":"IMPORTANT: Julia code is just-in-time compiled. This means that the first time you run a function, it will be compiled at that moment. So, please be patient! Future calls to the function will be much much faster. Trying out toy examples for the first calls is a good idea.","category":"page"},{"location":"man/installation/#Installation-of-the-PhyloClustering.jl-package","page":"Installation","title":"Installation of the PhyloClustering.jl package","text":"","category":"section"},{"location":"man/installation/","page":"Installation","title":"Installation","text":"To install the package, type inside Julia:","category":"page"},{"location":"man/installation/","page":"Installation","title":"Installation","text":"]\nadd PhyloClustering","category":"page"},{"location":"man/installation/","page":"Installation","title":"Installation","text":"The first step can take a few minutes, be patient.","category":"page"},{"location":"man/installation/","page":"Installation","title":"Installation","text":"The PhyloClustering.jl package has dependencies like Distributions and StatsBase (see the Project.toml file for the full list), but everything is installed automatically.","category":"page"},{"location":"man/installation/#Loading-the-Package","page":"Installation","title":"Loading the Package","text":"","category":"section"},{"location":"man/installation/","page":"Installation","title":"Installation","text":"To check that your installation worked, type this in Julia to load the package. This is something to type every time you start a Julia session:","category":"page"},{"location":"man/installation/","page":"Installation","title":"Installation","text":"using PhyloClustering","category":"page"},{"location":"man/installation/","page":"Installation","title":"Installation","text":"This step can also take a while, if Julia needs to pre-compile the code (after a package update for instance).","category":"page"},{"location":"man/installation/","page":"Installation","title":"Installation","text":"Press ? inside Julia to switch to help mode,  followed by the name of a function (or type) to get more details about it.","category":"page"}]
}
