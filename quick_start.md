# Quick start - for computerome users

This quick guide is for users of Computerome 2 - The Danish National Supercomputer for Life Sciences.

Go to the directory on computerome where you want to store your copy of flankophile.

```bash
git clone https://avthorn@bitbucket.org/genomicepidemiology/flankophile.git

cd flankophile

```

Prepare your 2 input files. A *Reference database* and a *Sample input list*. 
For details see [README](README.md). You can put these files in the input folder if you want.

Edit the [config file](config.yaml). Change the path of the input files so they refer to your input files. 
Change the values of other variables if you want to. See see [README](README.md) for more details.

Now use screen to to make a new session on computerome.

```bash
screen -S session_name

```

Now book and interactive node via the qsub system. Here is an example on how to do it.

```bash
qsub -W group_list=<<GROUPNAME>> -A <<GROUPNAME>> -l nodes=1:ppn=40,mem=20gb,walltime=14400 -I

```
Now you need to load the right modules:

```bash
module purge
module load tools  
module load miniconda3/4.11.0 
module load snakemake/6.9.1

```

Go to the flankophile folder. To run the pipeline:

```bash
snakemake --use-conda --cores 39

```

Now wait  for results.


**Problems with conda?**

If running the pipeline gives an error message related to conda then it 
may be because you do not have space for all the packages in your personal root folder.

As a default conda will download the packages needed to your root directory on Computerome. This root can only contain 10 GB and the packages for Flankophle takes up about 8 GB. You may need to delete some files in your personal root directory to make space. If you already have many packages in your hidden folder .conda, you may consider to delete it also `rm -r .conda`. 


After deleting some files then try to run the pipeline again.

If you want conda to store your packages and enviroments in another location than your root directoty then you can write this in the file called .condarc. It may be in your root directory. If you do not have the file .condarc then you can make one and conda will find it in your root.

```bash
auto_activate_base: false

pkgs_dirs:
    - /home/projects/SOMEPATH_YOU_CHOOSE/.conda/pkgs

envs_dirs:
    - /services/tools/miniconda3/4.11.0/envs
    - /home/projects/SOMEPATH_YOU_CHOOSE/.conda/envs

```



[Read the condarc_guide](https://conda.io/projects/conda/en/latest/user-guide/configuration/use-condarc.html#creating-and-editing)




**Visualisation with R studio**

R Studio on Computerome does **not** currently have the right version of R (June, 2022). You must run the R script on your local computer. Download the R folder to your computer. In the folder that contains the R folder you make a folder called output. In this folder you download folders 2_filter_gene_observations and 5_gene_clusters. Or all results foldes if you have infinite space on your computer. Now run the R script from the R folder. *Remember to use the correct version of R. Try 4.1.2 or 4.1.3.* Save plots as PDF. You can make the PDF longer than A4 if the results are too cramped. If you reach the maximum length that a pdf can be then save as png.
