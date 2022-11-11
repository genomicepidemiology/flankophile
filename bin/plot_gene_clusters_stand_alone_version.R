
# Companion script for Flankophile. 
# By Alex Vincent Thorn

# Run the plot in RStudio. R version 4.1.2 was used for development. Try with 4.1.2 or 4.1.3.  

# See all 4 plots by using the arrows in the plot viewer.
# The plots may look strange in the plot viewer. Export them to see how they really look.

# Save a plot as pdf in R studio from the plot viewer under "Export". Try A4 format. 
# If plot is ugly try to export as a pdf that is longer than A4. 

# Clusters with less than 3 gene observations will cause an error when run. Try another cluster.




### Settings ######################################################################################################

cluster_name <- "23_aph"   # Enter the name of a subfolder in output folder 4_cluster_results.

path_to_output_folder <- "../output/4_cluster_results/"   # relative path the the 4_cluster_results folder


# Choose desired legend position. "none" or "bottom". Use "None" if plot is very big, else use "bottom" to see the legend.
legendpos <- "bottom" 
#legendpos <- "none" 

x_lim_plot <- 0.7     # Adjust x-lim to make tiplaps visible in flank tree.

cutoff_limit_gene_names <- 11  # Try with different cutoffs.

############################################################################################################################

library(tidyverse)
library(ggtree)
library(gggenes)
library(treeio)
library(ape)
library(ggnewscale)


clus_gene <- str_split_fixed(cluster_name, pattern = "_", n=2)[,2]
clus_num <- str_split_fixed(cluster_name, pattern = "_", n=2)[,1]


cluster_results <- read_tsv(paste0(path_to_output_folder, cluster_name, "/", cluster_name, ".tsv"))



info_gene_ID <- cluster_results %>% 
  mutate(id = Gene_observation_ID) %>% 
  relocate(id) %>% 
  rename(IDENTITY = "%IDENTITY") %>% 
  mutate(Gene = str_c(GENE, IDENTITY, "%", sep = "_")) %>% 
  select(id, Gene) %>% 
  column_to_rownames(., var = "id") 


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


# Dealing with gene annotations ##################

cluster_prokka_raw <- read_tsv(paste0(path_to_output_folder, cluster_name, "/", cluster_name, ".gggenes"))


cluster_prokka <- cluster_prokka_raw %>% 
  mutate(gene_long = gene) %>% 
  mutate(gene = str_sub(gene, 1, cutoff_limit_gene_names))


# Find most common gene - the plot will be centered on this gene
most_commen_gene <- cluster_prokka %>% 
  group_by(gene) %>% 
  count(gene) %>% 
  arrange(desc(n))
most_commen_gene <- as.character(most_commen_gene[1,1])


# Manually force what gene to center on?
#most_commen_gene <- "my gene"

# Based on flank + gene sequence
target_and_flanking_regions.df = parseDistanceDF(paste0(path_to_output_folder, cluster_name, "/", cluster_name, ".target_and_flanking_regions_dist"))

# Based on flank sequence
flanks_masked_gene.df = parseDistanceDF(paste0(path_to_output_folder, cluster_name, "/", cluster_name, ".flanking_regions_only_dist"))

# Based on gene sequence
target_sequence_only.df = parseDistanceDF(paste0(path_to_output_folder, cluster_name, "/", cluster_name, ".target_sequence_only_dist"))


# Calculate tree from distance matrix
dist2tree = function(distmat) {
  clust = hclust(as.dist(distmat), method = "average")
  tree = as.phylo(clust)
  return(tree)
}


# Make trees

target_and_flanking_regions.tree = dist2tree(target_and_flanking_regions.df)
flanks_flanking_regions_only.tree = dist2tree(flanks_flanking_regions_only.df)
target_sequence_only.tree = dist2tree(target_sequence_only.df)


# flank + gene plot - based on gene plus flanking sequence ##################

p_target_and_flanking_regions <- ggtree(target_and_flanking_regions.tree, options(ignore.negative.edge=TRUE)) + 
  geom_tiplab(size = 1) +
  xlim_tree(x_lim_plot) +
  ggtitle(paste0("Cluster ", clus_num, " - ", clus_gene, " - distance tree based on gene plus flanking region sequences")) +
  geom_facet(mapping = aes(xmin = start, xmax = end, fill = gene, forward = orientation),
             data = cluster_prokka, geom = geom_motif, panel = 'Alignment',
             on = most_commen_gene, label = 'gene_long', align = 'left') +
  theme(legend.position=legendpos)

target_and_flanking_regions_plot <-facet_widths(p_target_and_flanking_regions, widths=c(1,2))

target_and_flanking_regions_plot



# flank  plot - based on only flanking sequence ##################

p_flanks_flanking_regions_only <- ggtree(flanks_flanking_regions_only.tree, options(ignore.negative.edge=TRUE)) + 
  geom_tiplab(size = 1) +
  xlim_tree(x_lim_plot) +
  ggtitle(paste0("Cluster ", clus_num, " - ", clus_gene, " - distance tree based on only flanking region sequences")) +
  geom_facet(mapping = aes(xmin = start, xmax = end, fill = gene, forward = orientation),
             data = cluster_prokka, geom = geom_motif, panel = 'Alignment',
             on = most_commen_gene, label = 'gene_long', align = 'left') +
  theme(legend.position=legendpos)

flanks_flanking_regions_only_plot <-facet_widths(p_flanks_flanking_regions_only, widths=c(1,2))

flanks_flanking_regions_only_plot


# gene  plot - based on only gene sequence ##################

p_target_sequence_only <- ggtree(target_sequence_only.tree, options(ignore.negative.edge=TRUE)) + 
  geom_tiplab(size = 1) +
  xlim_tree(x_lim_plot) +
  ggtitle(paste0("Cluster ", clus_num, " - ", clus_gene, " - distance tree based on only gene sequences")) +
  geom_facet(mapping = aes(xmin = start, xmax = end, fill = gene, forward = orientation),
             data = cluster_prokka, geom = geom_motif, panel = 'Alignment',
             on = most_commen_gene, label = 'gene_long', align = 'left') +
  theme(legend.position=legendpos)

target_sequence_only_plot <-facet_widths(p_target_sequence_only, widths=c(1,2))

target_sequence_only_plot


##Plot with reference gene % identity panel ###############################################################



# Change to target_sequence_only.tree or target_and_flanking_regions if wanted
tree <- flanks_flanking_regions_only.tree



t1 <- ggtree(tree, options(ignore.negative.edge=TRUE)) + geom_tiplab(size = 1)   


t2 <- facet_plot(t1, mapping = aes(xmin = start, xmax = end, fill = gene, forward = orientation),
                 data = cluster_prokka, geom = geom_motif, panel = 'Alignment',
                 on = most_commen_gene, label = 'gene_long', align = 'left') +
  labs(fill = "Gene") + 
  #scale_fill_brewer(name = "Gene", palette = "Set3") + new_scale_fill()
  scale_fill_hue() + new_scale_fill()



t3 <- gheatmap(t2, info_gene_ID,  width=0.1, hjust=0,
               colnames=TRUE, offset = 0, font.size = 1, colnames_position = "bottom", colnames_angle = 0) +
  labs(fill = "Target - Reference DB % identity") + 
  #scale_fill_brewer(name = "% Identity - Reference genes", palette = "Set1") + new_scale_fill()
  scale_fill_hue() + new_scale_fill()

t4 <- t3 + ggtitle(paste0("Cluster ", clus_num, " - ", clus_gene, " - distance tree based on only flanking region sequences")) + 
  new_scale_fill()

t4

#To save plots use ggsave

#ggsave(file = paste0("flanks_flanking_regions_only-ID-", cluster_name, "_tree", ".pdf"), plot = t4, width = 8.27, height = 11.69)
