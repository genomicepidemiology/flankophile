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
may be because you already have an enviroment loaded you need to delete. 

On Computerome you only have space for 10 GB in your personal root directory and the conda environment takes up some space.
Go to your personal root `cd ` and delete .conda: `rm -r .conda`. Also delete some other files if you have a lot of junk in your personal dir.
In the flankophile pipeline directory delete ./snakemake/conda and ./snakemake/conda-archieve.
Then try running the pipeline again.
