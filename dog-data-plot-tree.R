library(phytools)
trees <- read.tree("dog-data/dogs_estimated_gene_trees_merged6_cleaned.txt")
for (i in 1:length(trees)){
  trees[[i]] = reroot(trees[[i]], which(trees[[i]]$tip.label=="1"))
}
write.tree(tree, file='dogs_estimated_gene_trees_merged6_cleaned_rooted.txt') 


trees <- read.tree("dog-data/100k_taxa_cluster_1.trees")
for (i in 1:length(trees)){
  trees[[i]] = reroot(trees[[i]], which(trees[[i]]$tip.label=="Africanhuntingdog"))
}

writeNexus(trees, file="dog-data/100k_taxa_cluster_rooted_1.trees")

trees <- read.tree("dog-data/100k_taxa_cluster_2.trees")
for (i in 1:length(trees)){
  trees[[i]] = reroot(trees[[i]], which(trees[[i]]$tip.label=="Africanhuntingdog"))
}

writeNexus(trees, file="dog-data/100k_taxa_cluster_rooted_2.trees")

trees <- read.tree("dog-data/200k_taxa_cluster_1.trees")
for (i in 1:length(trees)){
  trees[[i]] = reroot(trees[[i]], which(trees[[i]]$tip.label=="Africanhuntingdog"))
}

writeNexus(trees, file="dog-data/200k_taxa_cluster_rooted_1.trees")

trees <- read.tree("dog-data/200k_taxa_cluster_2.trees")
for (i in 1:length(trees)){
  trees[[i]] = reroot(trees[[i]], which(trees[[i]]$tip.label=="Africanhuntingdog"))
}

writeNexus(trees, file="dog-data/200k_taxa_cluster_rooted_2.trees")