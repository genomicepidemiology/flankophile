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

Now use screen to to make a new session on Computerome.

```bash
screen -S session_name

```

Now book and interactive node via the qsub system. Here is an example on how to do it.

```bash
qsub -W group_list=<<GROUPNAME>> -A <<GROUPNAME>> -l nodes=1:thinnode:ppn=40,mem=180gb,walltime=25200 -I
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


**Problems with conda?**

If running the pipeline gives an error message related to conda and packages then it 
may be because you do not have space for all the packages in your personal root folder.

As a default conda will download the packages needed to your root directory on Computerome. This root can only contain 10 GB and the packages for Flankophle takes up more than 8 GB. 

You may need to force  conda to store your packages and enviroments in another location than your root directoty. You can write where in the file called ".condarc". It may be in your root directory. If you do not have the file ".condarc" in then you can make one and place it in your root and conda will find it.

Here is what you need to write:

```bash
auto_activate_base: false

pkgs_dirs:
    - /home/projects/SOMEPATH_YOU_CHOOSE/.conda/pkgs

envs_dirs:
    - /home/projects/SOMEPATH_YOU_CHOOSE/.conda/envs

```



[Read the condarc_guide](https://conda.io/projects/conda/en/latest/user-guide/configuration/use-condarc.html#creating-and-editing)





