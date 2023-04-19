

# R script for Flankophile. 
# By Alix Vincent Thorn

# R version 4.1.2 or 4.1.3.  


library(tidyverse)
library(ggtree)       # ggtree_3.2.1
library(gggenes)      # gggenes_0.4.1
library(treeio)       # treeio_1.18.1
library(ape)          # ape_5.6-2
library(ggnewscale)   # ggnewscale_0.4.7


##########################################################################################################


path_to_output_folder <- "output"   # relative path to the entire flankophile output folder


path_to_folder_for_plots <- paste0(path_to_output_folder, "/4_plots") 



path_to_all_cluster_folders <- paste0(path_to_output_folder, "/4_cluster_results/")

path_to_config_file <- paste0("config.yaml")

config <- read_tsv(path_to_config_file, col_names = FALSE)

upstream <- config %>% filter(str_starts(X1, "flank_length_upstreams:", negate = FALSE)) %>% 
  mutate(X1 = str_extract(X1, "[:digit:]+")) 

upstream <- as.integer(upstream[1,1])




##### #####


make_plots <- function(cluster_name) {
  
  
    
  cluster_results <- read_tsv(paste0(path_to_all_cluster_folders, cluster_name, "/", cluster_name, ".tsv"), col_types = cols(GENE = col_character(),
                                                                                                                           METADATA = col_character())) %>% 
    replace_na(list(METADATA = "NA"))
 
  cluster_prokka_raw <- read_tsv(paste0(path_to_all_cluster_folders, cluster_name, "/", cluster_name, ".gggenes"), col_types = cols(molecule = col_character(),
                                                                                                                                    gene = col_character(),
                                                                                                                                    start = col_integer(),
                                                                                                                                    end = col_integer(),
                                                                                                                                    orientation = col_integer(),
                                                                                                                                    db_name = col_character()))
  
  
  if (length(unique(cluster_results$METADATA)) == 1 & unique(cluster_results$METADATA)[1] == "No_metadata") {
  metadata_present = "no"
  } else {
  metadata_present = "yes"
  }
  

  
  
  length_of_target <- cluster_results %>% mutate(len_gen = END - START + 1) %>%
    group_by(len_gen) %>% 
    count(len_gen) %>% 
    arrange(desc(n))
  
  length_of_target <- as.double(length_of_target[1,1])
  
  midpoint <- length_of_target / 2 + upstream 
  
  
    # Dealing with gene annotation names ##################
  
  cutoff_limit_gene_arrow_labels <- 8  # Try with shorter value if some genes labels are not shown on the arrow.
  
  
  cluster_prokka <- cluster_prokka_raw %>%
    mutate(Gene = case_when(start < midpoint & end > midpoint ~ "TARGET",
                            TRUE ~ gene)) %>% 
    mutate(gene_arrow_label = case_when(gene == "Hypothetical protein" ~ "",
                                        nchar(gene) < cutoff_limit_gene_arrow_labels ~ gene,
                                        nchar(gene) >= cutoff_limit_gene_arrow_labels ~ paste0(str_sub(gene, 1, cutoff_limit_gene_arrow_labels), "...")))
    
  # Deal with gene arrow colours
  
  gg_color_hue <- function(n) {
    hues = seq(30, 330, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  
  
  target_hp_col_vec <- c("#FF0000", "#D3D3D3")
  names(target_hp_col_vec) <- c("TARGET", "Hypothetical protein")
  

  novelty_gene_table <- cluster_prokka %>% select(Gene) %>%  filter(Gene != "TARGET" & Gene != "Hypothetical protein")         
  
  novelty_genes <- unique(novelty_gene_table$Gene)
  
  novelty_genes_vec <- gg_color_hue(length(novelty_genes))
  
  names(novelty_genes_vec) <- novelty_genes
  
  
  all_genes_vec <- append(target_hp_col_vec, novelty_genes_vec)
  
  #scale_fill_manual(values = all_genes_vec)
  
  
  # Deal with a completely empty plot ;(
    
  if (nrow(cluster_prokka) == 0) {    
    t4 <- ggplot() +
      theme_void() +
      geom_text(aes(0,0,label='No genes annotated by Prokka')) +
      xlab(NULL) +
      ggtitle(paste0("Cluster ", cluster_name))
    
    t4
    ggsave(file = paste0(path_to_folder_for_plots, "/", cluster_name, ".pdf"), plot = t4, width = 210, height = 150, units = "mm")
    
  } 
  
  
  
  if (nrow(cluster_results) == 1 & nrow(cluster_prokka) > 0) {    
    t4 <- ggplot(cluster_prokka, aes(xmin = start, xmax = end, y = molecule, fill = Gene, forward = orientation)) +
      geom_gene_arrow() +
      facet_wrap(~ molecule, scales = "free", ncol = 1) +
      #scale_fill_brewer(palette = "Set3") +
      #scale_fill_hue(direction = 1) + 
      scale_fill_manual(values = all_genes_vec)
      ggtitle(paste0("Cluster ", cluster_name))
    
    t4
    ggsave(file = paste0(path_to_folder_for_plots, "/", cluster_name, ".pdf"), plot = t4, width = 210, height = 150, units = "mm")
    
  } 
  
  
  
  
    
  info_meta <- cluster_results %>% 
  mutate(id = OBSERVATION_ID) %>% 
  select(id, METADATA) %>%
  column_to_rownames(., var = "id") 
    
  
  
  
  
  info_gene_ID <- cluster_results %>% 
    mutate(id = OBSERVATION_ID) %>% 
    relocate(id) %>% 
    rename(IDENTITY = "%IDENTITY") %>% 
    mutate(vari = str_c(GENE, IDENTITY, "%", VARIANT, sep = "_")) %>% 
    select(id, vari) %>%
    rename(VARIANT = vari) %>%
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
  
  

  
  
  
  
  # Based on flank + gene sequence
  target_and_flanking_regions.df = parseDistanceDF(paste0(path_to_all_cluster_folders, cluster_name, "/", cluster_name, ".target_and_flanking_regions_dist"))
  
  # Based on flank sequence
  flanks_flanking_regions_only.df = parseDistanceDF(paste0(path_to_all_cluster_folders, cluster_name, "/", cluster_name, ".flanking_regions_only_dist"))
  
  # Based on gene sequence
  target_sequence_only.df = parseDistanceDF(paste0(path_to_all_cluster_folders, cluster_name, "/", cluster_name, ".target_sequence_only_dist"))
  
  
  # Calculate tree from distance matrix
  dist2tree = function(distmat) {
    clust = hclust(as.dist(distmat), method = "average")
    tree = as.phylo(clust)
    return(tree)
  }
  
  
  if (nrow(cluster_results) >= 2) {
    # Make trees
    
    target_and_flanking_regions.tree = dist2tree(target_and_flanking_regions.df)
    flanks_flanking_regions_only.tree = dist2tree(flanks_flanking_regions_only.df)
    target_sequence_only.tree = dist2tree(target_sequence_only.df)
    
    ##Plot  ###############################################################
    plot_height <- 150 + (nrow(cluster_results) * 4) + (length(unique(cluster_results$METADATA)) * 10)
    

    
    #Flank only
    t1 <- ggtree(flanks_flanking_regions_only.tree, options(ignore.negative.edge=TRUE)) + geom_tiplab(size = 2)
    
    off <- max(t1$data$x) / 9
    
    off2 <- off * 2
    
    
    
    t2 <- facet_plot(t1, mapping = aes(xmin = start, xmax = end, fill = Gene, forward = orientation),
                     data = cluster_prokka, geom = geom_motif, panel = 'Alignment',
                     on = "TARGET", label = 'gene_arrow_label', align = 'left') +
      labs(fill = "Gene") + 
      scale_fill_manual(values = all_genes_vec) +
      #scale_fill_hue(direction = 1) + 
      new_scale_fill()
    
    
    t3 <- gheatmap(t2, info_gene_ID,  width=0.1, hjust=1,
                   colnames=TRUE, offset = off, font.size = 1, colnames_position = "bottom", colnames_angle = 0) +
      labs(fill = "Target variant") + 
      scale_fill_hue(direction = -1, l = 70, c = 30) + new_scale_fill() +
      theme(legend.text = element_text(size=5))
      
      
    if (metadata_present == "yes") {
      tx <- gheatmap(t3, info_meta,  width=0.1, hjust=1,
                   colnames=TRUE, offset = off2 , font.size = 1, colnames_position = "bottom", colnames_angle = 0) +
        labs(fill = "Metadata") + 
         scale_fill_hue(direction = -1, l = 85) + new_scale_fill() +
        theme(legend.text = element_text(size=5))
      t4 <- tx + ggtitle(paste0("Cluster ", cluster_name, " - distance tree based on flanking region sequences only")) + 
        new_scale_fill()  + if(length(unique(info_gene_ID$VARIANT)) > 25){ theme(legend.position="none")}
    }  else {
      t4 <- t3 + ggtitle(paste0("Cluster ", cluster_name, " - distance tree based on flanking region sequences only")) + 
        new_scale_fill() + if(length(unique(info_gene_ID$VARIANT)) > 25){ theme(legend.position="none")}
    }
    
    
    t4
    
    ggsave(file = paste0(path_to_folder_for_plots, "/", cluster_name, "_", "flanking_regions_only", ".pdf"), plot = t4, width = 210, height = plot_height, units = "mm", limitsize = FALSE)
    
    
    
    ## Gene only
    
    t1 <- ggtree(target_sequence_only.tree, options(ignore.negative.edge=TRUE)) + geom_tiplab(size = 2)   
    
    off <- max(t1$data$x) / 9
    
    off2 <- off * 2
    
    
    t2 <- facet_plot(t1, mapping = aes(xmin = start, xmax = end, fill = Gene, forward = orientation),
                     data = cluster_prokka, geom = geom_motif, panel = 'Alignment',
                     on = "TARGET", label = 'gene_arrow_label', align = 'left') +
      labs(fill = "Gene") + 
      scale_fill_manual(values = all_genes_vec) +
      #scale_fill_hue(direction = 1) + 
      new_scale_fill()
    
    
    t3 <- gheatmap(t2, info_gene_ID,  width=0.1, hjust=1,
                   colnames=TRUE, offset = off, font.size = 1, colnames_position = "bottom", colnames_angle = 0) +
      labs(fill = "Target variant") + 
      scale_fill_hue(direction = -1, l = 70, c = 30) + new_scale_fill() +
      theme(legend.text = element_text(size=5))
      
      
    if (metadata_present == "yes") {
      tx <- gheatmap(t3, info_meta,  width=0.1, hjust=1,
                   colnames=TRUE, offset = off2 , font.size = 1, colnames_position = "bottom", colnames_angle = 0) +
        labs(fill = "Metadata") + 
         scale_fill_hue(direction = -1, l = 85) + new_scale_fill() +
        theme(legend.text = element_text(size=5))
      t4 <- tx + ggtitle(paste0("Cluster ", cluster_name, " - distance tree based on target sequences only")) + 
        new_scale_fill()  + if(length(unique(info_gene_ID$VARIANT)) > 25){ theme(legend.position="none")}
    }  else {
      t4 <- t3 + ggtitle(paste0("Cluster ", cluster_name, " - distance tree based on target sequences only")) + 
        new_scale_fill() + if(length(unique(info_gene_ID$VARIANT)) > 25){ theme(legend.position="none")}
    }  
    
 
    t4
    
    ggsave(file = paste0(path_to_folder_for_plots, "/", cluster_name, "_", "target_sequence_only", ".pdf"), plot = t4, width = 210, height = plot_height, units = "mm", limitsize = FALSE)
    
    
    # Both flank and gene
    
    t1 <- ggtree(target_and_flanking_regions.tree, options(ignore.negative.edge=TRUE)) + geom_tiplab(size = 2)   
    
    off <- max(t1$data$x) / 9
    
    off2 <- off * 2
    
    
    t2 <- facet_plot(t1, mapping = aes(xmin = start, xmax = end, fill = Gene, forward = orientation),
                     data = cluster_prokka, geom = geom_motif, panel = 'Alignment',
                     on = "TARGET", label = 'gene_arrow_label', align = 'left') +
      labs(fill = "Gene") + 
      scale_fill_manual(values = all_genes_vec) +
      #scale_fill_hue(direction = 1) + 
      new_scale_fill()
    
    
    t3 <- gheatmap(t2, info_gene_ID,  width=0.1, hjust=1,
                   colnames=TRUE, offset = off, font.size = 1, colnames_position = "bottom", colnames_angle = 0) +
      labs(fill = "Target variant") + 
      scale_fill_hue(direction = -1, l = 70, c = 30) + new_scale_fill() +
      theme(legend.text = element_text(size=5))
    
    
    if (metadata_present == "yes") {
      tx <- gheatmap(t3, info_meta,  width=0.1, hjust=1,
                   colnames=TRUE, offset = off2 , font.size = 1, colnames_position = "bottom", colnames_angle = 0) +
        labs(fill = "Metadata") + 
         scale_fill_hue(direction = -1, l = 85) + new_scale_fill() +
        theme(legend.text = element_text(size=5))
      t4 <- tx + ggtitle(paste0("Cluster ", cluster_name, " - distance tree based on target and flanking region sequences")) + 
        new_scale_fill()  + if(length(unique(info_gene_ID$VARIANT)) > 25){ theme(legend.position="none")}
    }  else {
      t4 <- t3 + ggtitle(paste0("Cluster ", cluster_name, " - distance tree based on target and flanking region sequences")) + 
        new_scale_fill() + if(length(unique(info_gene_ID$VARIANT)) > 25){ theme(legend.position="none")}
    } 
    
    
    
    t4
    
    ggsave(file = paste0(path_to_folder_for_plots, "/", cluster_name, "_", "target_and_flanking_regions", ".pdf"), plot = t4, width = 210, height = plot_height, units = "mm", limitsize = FALSE)
  }
}


all_cluster_names <- dir(path = path_to_all_cluster_folders, 
                         pattern = "*") 

for(cluster in all_cluster_names) {
  
  make_plots(cluster)
}


# Run just one cluster
#make_plots("99_catB2_1_AF047479")   # Write cluster name






