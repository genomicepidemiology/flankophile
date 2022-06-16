# FLANKOPHILE
FLANKOPHILE version 0.0.4
By Alex Thorn

## About Flankophile

Flankophile is a pipeline build for easy analysis of gene synthony - the genetic context of genes. Flankophile locates genes of interest in the input data and extracts the gene plus the upstream and downstream flanking regions.

Requied input: Assempled genetic data from multiple samples in fasta format plus a reference gene database in fasta format.

Output: Distance matrix, neighbor joining trees, genetic sequences of genes and flanking regions. 

Flankophile is a Snakemake based pipeline. Snakemake is a python based workflow management system.
The [Snakefile](Snakefile) contains the main pipeline code.

Flankofile is conda based and uses a conda env. You can see it in the file environment.yaml.
When you run the pipeline Snakemake and conda will download the required conda tools.  

## How to run Flankophile

Download this repository to your computer.

```bash
git clone https://avthorn@bitbucket.org/genomicepidemiology/flankophile.git

```


### Prerequisites


**If you are using Computerome 2**

If you are using Computerome 2 (The Danish National Supercomputer for Life Sciences) then you can use the module system. 
Book an interactive node via the qsub system. 
Load these modules:  
 
`module load tools`  
`module load miniconda3/4.11.0`  
`module load snakemake/6.9.1`


**If you are not using Computerome 2**

You need to miniconda and Snakemake to run Flankophile.



### Input files

In order to use the pipeline you need to prepare two input files: The reference database and the data you want to analyse.
 
The path to the input files must be given in the config file [config.yaml](config.yaml).


#### Reference database

The reference database contains reference sequences of all the genes that you want to search for in your data. 
Flankophile will find these genes in your data and extract the gene sequence and the flanking sequences on each side of the gene.
The reference database can contain reference genes that are not homologs. Flankophile will cluster the reference genes by percentage identity and report the results separately for each cluster.
The reference database must be a multifasta file. The resfinder database is given as an example of a reference database.

It can be found here [Resfinder_08_02_2022_dub_rem.fa](input/example_input_files/Resfinder_08_02_2022_dub_rem.fa). The ResFinder database consists 
of acquired antimicrobial resistance genes. The version found in this repository is from February 8 2022. 
The uptodate ResFinder database is found [here](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/).


#### Sample input list

Your sample input data must consist of a number of assemblies or genemes in multifasta format. One multifasta per sample. You can input as many samples as wanted. You have to use either an assembly or a contig input file.

**Assembly level analysis**

The recommended way is to use assembly-level analysis. This is when you want to analyse the entire assembly and not analyse only a subset of the contigs in each assembly. An assembly can also be an individual bin from binned data, or it can be a fasta of unbinned contigs. The input file is a tsv file with two columns. The first column has to be a unique name for each input fasta, for example, "sample_1" or "e.coli_bin_32". 

The second column is the full path to the fasta file, including the file name. The pipeline will ignore lines that start with #. This is useful if you want to add human-readable headers between different datasets in the file.
In [test_input_assemblies_list_small.tsv](input/example_input_files/test_input_assemblies_list_small.tsv) you can see an example of an assembly level input file. 

**Config level analysis**

Three columns. One line per contig. The first column has to be a unique name for each input fasta. The second is the name of the contig. The third column is the full path to the fasta. 
In [test_input_contig_list_small.tsv](input/example_input_files/test_input_contig_list_small.tsv) you can see an example of a contig level input file. 


### Configuration file

Fill out [config.yaml](config.yaml) before running the pipeline. The path to input files can be given as relative path to the location of the Snakefile.
The configuation file contains numbered sections. Each number refer to an output folder.


| **Variable name**                   | **Suggestion**        | **Variable**    | **Notes**                                                                                          |
|-------------------------------------|-----------------------|-----------------|----------------------------------------------------------------------------------------------------|
| database                            | "input/db.fa"         | Path to file.   | d                                                                                                  |
| input_format                        | "assemblies"          | Format.         | d                                                                                                  |
| input_list                          | "input/sa.tsv"        | Path to file.   | d                                                                                                  |
| flank_length_upstreams              | "3000"                | Length in bp.   | d                                                                                                  |
| flank_length_downstreams            | "3000"                | Length in bp.   | d                                                                                                  |
| min_coverage_abricate               | "95"                  |                 | d                                                                                                  |
| min_identity_abricate               | "95"                  | In percent.     | d                                                                                                  |
| cluster_identity_cd_hit             | "0.95"                |                 | See https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#CDHITEST                              |
| cluster_wordsize_cd_hit             | "9"                   |                 | See https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#CDHITEST                              |
| cluster_length_dif_cd_hit           | "0.9"                 |                 | See https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#CDHITEST                              |
| Kmersize_kma                        | "16"                  | Kmer size.      | For kma index.                                                                                     |
| distance_measure                    | "1"                   | Distmatrix      | Methods: 1 k-mer hamming distance 64 Jaccard distance 256 Cosine distance 4096 Chi-square distance |


### Running the pipeline

Enter the flankophile folder.

Run the pipeline: 
`snakemake --use-conda --cores 39` 

Cores is the number of cores available. For more info on flags visit: 
https://snakemake.readthedocs.io/en/stable/executing/cli.html#command-line-interface To add more samples. Delete folders and files numbered 2 and up.


#### Running only part of the pipeline

If your dataset is very large you may want to run the pipeline in 3 stages. That way you may avoid
spending compute time rerunning the later stages again and again while you parametrize the first stages.
In the Snakefile around line 30 you find *rule all*. There are 4 files under the rule. 
The first and last file must always be included. 

If you write a # in front of "output/99_trees_and_distance_matrixes_done" 
then you do not run the last stage, step 5 of the pipeline. You can remove the # again later to run this step.
If you write a # in front of both
"output/4_cluster_by_gene_family/cd_hit.processed" and "output/99_trees_and_distance_matrixes_done"
 then you will only run the first stage, step 1 and 2 of the pipeline.

#### Rerunning the pipeline
If you want to rerun part of the pipeline you simply delete output folders and files and 
Snakemake will rerun these part of the pipeline when you give the run command `snakemake --use-conda --cores 39` again. 
You must delete output in descending order starting with 99 and down to where you do not want to rerun anymore. 

**New config settings**

If you want to rerun the pipeline with new settings you must look in the config file [config.yaml](config.yaml)
 and take note of the numbered sections. If you want to change a parameters you must delete all output
 with a number equal or higher than the numbe of the section. 

**More data**

If you want to add more data you must delete all output except 
0_setup_abricate_db and 1_abricate. You add more data by appending more lines to the
 input_list file. If you want to use a difference reference database you must delete the entire
 output directory and start over. 

#### Using Computerome 2 and having problems with conda?
If running the pipeline gives an error message related to conda then it 
may be because you already have an enviroment loaded you need to delete. 

On C2 you only have space for 10 GB in your personal root directory and the conda environment takes up some space.
Go to your personal root `cd ` and delete .conda: `rm -r .conda`. Also delete some other files if you have a lot of junk in your personal dir.
In the Flankophile pipeline directory delete ./snakemake/conda and ./snakemake/conda-archieve.



## Output

The output directories are created by Flankophile in the order of the numbers. 


**1_abricate**

Contains one directory per input fasta. Can be deleted if one is done with the analysis.

**2_filter_gene_observations**

In this directory abricate results from all the samples have been merged. It contains 3 reports and 3 tsv files.
In step 1 no filtering has taken place yet. The raw collective abricate results can be found in 1_abricate_all.tsv.
In step 2 the results have been filtered on minimum coverage and minimum sequence identity. 
2_abricate_filter_qual.report contains information on how many gene observtions that were discarded during filtering.

In step 3 the tsv from step 2 has been filtered on flank length. Gene observations that are so close to the
 edge of a contig that there is not space for the full desired flank length are discarded. 
 3_abricate_filter_length.report contains information on how many gene observtions that were discarded during filtering.
 3_final_gene_results.tsv contains all the gene observations that are included in the further analysis.


**3_extract_sequences**

These files are just temperaty files and can be deleted by the user without loss of important data.

**4_cluster_by_gene_family**

cd_hit.clstr contains the results of clustering the reference genes used as templates by identity.
The rest of the contents of this folder can be deleted by the user without loss of important data.
  
**5_gene_clusters**

Contains one directory for each reference gene cluster. Directory names have two parts. The first part is a unique number. 
The second part after the underscore is the first part of the name of the gene that seeded the cluster.  

## Visualization

The R script [plot_gene_clusters_from_flankophile.R](R/plot_gene_clusters_from_flankophile.R)
 can be used to visualize the results in output folder 5.
 It plots distance trees with gene annotation for a gene cluster. 
The trees are based on the distance matrices. 
The flank tree is based on the flanking regions alone since the gene is masked. 
The distance matrix file ends with .masked_gene_tree.
 
The gene tree is based on the target genes sequence alone. The distance matrix ends with .just_gene_dist. 
The script can be run with RStudio and should be run with the R directory as the working directory.

## Contact
Alex Vincent Thorn
alvit@food.dtu.dk


## License

## Acknowledgement


## Citations

This section is not yet complete!


**Abricate**

By Torsten Seemann

[https://github.com/tseemann/abricate](https://github.com/tseemann/abricate)




**any2fasta**


**bedtools**


**blast**


**cd-hit**


**clearcut**

[Clearcut: a fast implementation of relaxed neighbor joining](https://academic.oup.com/bioinformatics/article/22/22/2823/197536)



**clstr2txt.pl**

[https://raw.githubusercontent.com/weizhongli/cdhit/master/clstr2txt.pl](https://raw.githubusercontent.com/weizhongli/cdhit/master/clstr2txt.pl)




**kma**

By Philip T.L.C. Clausen
  
[https://bitbucket.org/genomicepidemiology/kma/src/master/](https://bitbucket.org/genomicepidemiology/kma/src/master/)

[Rapid and precise alignment of raw reads against redundant databases with KMA](https://pubmed.ncbi.nlm.nih.gov/30157759/)




**prokka**

By Torsten Seemann

[https://github.com/tseemann/prokka](https://github.com/tseemann/prokka)

[Prokka: rapid prokaryotic genome annotation](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517)




**seqkit**


**Snakemake**





