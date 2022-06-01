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

The output folders are created by Flankophile in the order of the numbers. 

### 2_filter_gene_observations

### 5_gene_clusters

### Visualization

## Citations

## License

## Contact
Alex Thorn
alvit@food.dtu.dk
