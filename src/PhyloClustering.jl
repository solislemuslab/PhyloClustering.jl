__precompile__()

module PhyloClustering
    include("bipartition.jl")
    include("helper.jl")
    include("models.jl")

    export 
        kmeans_label, 
        gmm_label, 
        hc_label, 
        dbscan_label,
        standardize_tree,
        plot_clusters,
        distance,
        num_bipartitions,
        show_bipartitions,
        show_bipartition,
        get_bipartition,
        print_bipartition,
        num_to_name
end