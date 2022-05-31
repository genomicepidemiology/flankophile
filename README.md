# FLANKOPHILE #
FLANKOPHILE version 0.0.3
By Alex Thorn


Flankofile is conda based and uses a conda env. You can see it in the file environment.yaml. Edit the config file config.yaml to your needs. You need to prepare an input file, a tsv. Read about it 
in the config file.


### Using FLANKOPHILE on Computerome ###
Use an interactive node. 
Load these modules: 
`module load tools` 
`module load miniconda3/4.11.0` 
`module load snakemake/6.9.1`

### Running the pipeline ###

You need miniconda to run the pipeline.

Git clone the whole repo to your computer.

Run the pipeline: 
`snakemake --use-conda --cores 39` 

Cores is the number of cores available. For more info on flags visit: 
https://snakemake.readthedocs.io/en/stable/executing/cli.html#command-line-interface To add more samples. Delete folders and files numbered 2 and up.

Problems with conda, may you already have an enviroment loaded you need to delete:
Delete .conda in root directory.
Delele ./snakemake/conda and ./snakemake/conda-archieve in snakemake directory.


### Contact ###
Alex Thorn
alvit@food.dtu.dk
