# FLANKOPHILE
FLANKOPHILE version 0.0.3
By Alex Thorn

## What does Flankophile do?


## How to run Flankophile

Flankophile is a Snakemake based pipeline. Snakemake is a python based workflow management system.
The [Snakefile](Snakefile) contains the main pipeline code.

### Prerequisites

You need miniconda to run the pipeline.

Flankofile is conda based and uses a conda env. You can see it in the file environment.yaml.
When you run the pipeline Snakemake and conda will download the required conda tools.  


#### If you are using Computerome 2 
If you are using Computerome 2 (The Danish National Supercomputer for Life Sciences) then you can use the module system. 
Book an interactive node via the qsub system. 
Load these modules:  
 
`module load tools`  
`module load miniconda3/4.11.0`  
`module load snakemake/6.9.1`

### Download Flankophile

### Input files


### Configuration file


### Running the pipeline


Run the pipeline: 
`snakemake --use-conda --cores 39` 

Cores is the number of cores available. For more info on flags visit: 
https://snakemake.readthedocs.io/en/stable/executing/cli.html#command-line-interface To add more samples. Delete folders and files numbered 2 and up.


#### Problems with conda?
If running the pipeline gives an error message related to conda then it 
may be because you already have an enviroment loaded you need to delete:
Delete .conda in your personal root directory. 
Delele ./snakemake/conda and ./snakemake/conda-archieve in snakemake directory.

#### Rerunning the pipeline

##### New config settings

##### More data

## Output

The output directories are created by Flankophile in the order of the numbers. 


### 1_abricate
Contains one directory per input fasta. Can be deleted if one is done with the analysis.

### 2_filter_gene_observations
In this directory abricate results from all the samples have been merged. It contains 3 reports and 3 tsv files.
In step 1 no filtering has taken place yet. The raw collective abricate results can be found in **1_abricate_all.tsv**.
In step 2 the results have been filtered on minimum coverage and minimum sequence identity. 
**2_abricate_filter_qual.report** contains information on how many gene observtions that were discarded during filtering.
In step 3 the tsv from step 2 has been filtered on flank length. Gene observations that are so close to the
 edge of a contig that there is not space for the full desired flank length are discarded. 
 **3_abricate_filter_length.report** contains information on how many gene observtions that were discarded during filtering.
 **3_final_gene_results.tsv** contains all the gene observations that are included in the further analysis.


### 3_extract_sequences
These files are just temperaty files and can be deleted by the user without loss of important data.

### 4_cluster_by_gene_family
**cd_hit.clstr** contains the results of clustering the reference genes used as templates by identity.
The rest of the contents of this folder can be deleted by the user without loss of important data.
  
### 5_gene_clusters
Contains one directory for each reference gene cluster. Directory names have two parts. The first part is a unique number. 
The second part after the underscore is the first part of the name of the gene that seeded the cluster.  

### Visualization

The R script [plot_gene_clusters_from_flankophile.R](R/plot_gene_clusters_from_flankophile.R) can be used to visualize the results
in output folder 5. It plots distance trees with gene annotation for a gene cluster. The script can be run with RStudio.

## Citations

## License

## Contact
Alex Thorn
alvit@food.dtu.dk
