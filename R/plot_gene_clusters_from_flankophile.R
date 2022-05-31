
# Companion script for Flankophile. Run script in RStudio from R folder.
# Save plot as pdf in R studio from the plot viewer. Make plot bigger than A4 during export if you cannot read the text. 

cluster_name <- "1_blaTEM-57"   # Enter the name of a subfolder in output folder 5_gene_clusters.

legendpos <- "bottom" # Enter desired legend position. "None" or "bottom". Use "None" if plot is very big.

x_lim_plot <- 0.7     # Adjust x-lim to make tiplaps visible in flank tree.

############################################################################################################################

library(tidyverse)
library(ggtree)
library(gggenes)
library(treeio)
library(ape)


path_to_kma_dist_matrix__masked_gene <- paste0("../output/5_gene_clusters/", cluster_name, "/", cluster_name, ".masked_gene_dist")

path_to_kma_dist_matrix__just_gene <- paste0("../output/5_gene_clusters/", cluster_name, "/", cluster_name, ".just_gene_dist")


# Function that imports distance matrix from KMA dist

parseDistanceDF = function(phylip_file) {
  # Read the first line of the phylip file to find out how many sequences/samples it contains
  temp_connection = file(phylip_file, 'r')
  len = readLines(temp_connection, n=1)
  len = as.numeric(len)
  len = len + 1
  close(temp_connection)
  # Make distance data frame
  phylip_data = read.delim(phylip_file, 
                           fill=T,
                           skip=0, 
                           col.names=1:len)
  seqNames = make.unique(phylip_data[,1])
  # Deal with odd sample names containing something after a space
  seqNames = str_split(seqNames, pattern = " ", simplify = T)[,1]
  rownames(phylip_data) = seqNames
  phylip_data = phylip_data[,-1]
  colnames(phylip_data) = row.names(phylip_data)
  return(phylip_data)
}




# Based on flank sequence
flank.df = parseDistanceDF(path_to_kma_dist_matrix__masked_gene)


# Based on gene sequence
gene.df = parseDistanceDF(path_to_kma_dist_matrix__just_gene)


# Calculate tree from distance matrix
dist2tree = function(distmat) {
  clust = hclust(as.dist(distmat), method = "average")
  tree = as.phylo(clust)
  return(tree)
}


# Make trees

flank.tree = dist2tree(flank.df)
gene.tree = dist2tree(gene.df)



# gene plot ####################
gene_plot <- ggtree(gene.tree) + 
  geom_tiplab() +
  ggtitle(paste0(cluster_name, " - gene"))

gene_plot

# flank plot ##################

cluster_prokka_raw <- read_tsv(paste0("../output/5_gene_clusters/", cluster_name, "/", cluster_name, ".gggenes"))

cluster_prokka <- cluster_prokka_raw %>% 
  rename( gene_name_full = gene) %>% 
  mutate(gene = case_when(str_count(gene_name_full, "_") >= 1 ~ str_split(gene_name_full, pattern = "_", simplify = T)[,1],
                          str_count(gene_name_full, "_") < 1 ~ gene_name_full))

# Find most common gene
most_commen_gene <- cluster_prokka %>% 
  group_by(gene) %>% 
  count(gene) %>% 
  arrange(desc(n))
most_commen_gene <- as.character(most_commen_gene[1,1])

p <- ggtree(flank.tree) + 
  geom_tiplab() +
  xlim_tree(x_lim_plot) +
  ggtitle(paste0(cluster_name, " - flanking regions")) +
  geom_facet(mapping = aes(xmin = start, xmax = end, fill = gene, forward = orientation),
             data = cluster_prokka, geom = geom_motif, panel = 'Alignment',
             on = most_commen_gene, label = 'gene', align = 'left') +
  theme(legend.position=legendpos)

flank_plot <-facet_widths(p, widths=c(1,2))

flank_plot



# R version 4.1.2
