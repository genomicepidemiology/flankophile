# By Alex Vincent Thorn


import os

import sys

import argparse

import re

gff_path = str(sys.argv[1])

gggenes_path = str(sys.argv[2])

abricate_path = str(sys.argv[3])


config_path = "config.yaml"


gff_file = open(gff_path,"r")

gggenes_file = open(gggenes_path,"w")




config_file = open(config_path,"r")

for con_line in config_file:
    if con_line.startswith("flank_length_upstreams:"):
        uplength = int(con_line.split('"')[1])

config_file.close()


outheader = "molecule" + "\t" + "gene" + "\t" + "start" + "\t" + "end" + "\t" + "orientation" + "\t" + "db_name"

print(outheader, file=gggenes_file)


abricate_dict = {}

abricate_file = open(abricate_path,"r")

for a_line in abricate_file:
    if not a_line.startswith("#"):
        
        a_line_list = a_line.split("\t")
        
        a_molecule = a_line_list[14]
        
        a_gene = a_line_list[5]
        
        a_start = uplength + 1
        
        a_end = int(a_line_list[3]) - int(a_line_list[2]) + 1 + uplength
             
        a_orientation = 1   #always forward since we made the fasta file so
               
        
        a_db_name = "user_db"
        
        abricate_dict[a_molecule] = [a_start, a_end]
    
        a_outline = str(a_molecule) + "\t" + str(a_gene) + "\t" + str(a_start) + "\t" + str(a_end) + "\t" + str(a_orientation) + "\t" + str(a_db_name)
        print(a_outline, file=gggenes_file)



fasta_flag = False

for line in gff_file:
    line = line.strip()
    if line.startswith("##FASTA"):
        fasta_flag = True
    if fasta_flag == False and not line.startswith("##"):
        line_list = line.split()
        molecule = line_list[0]
        start = int(line_list[3])
        end = int(line_list[4])
        strand = line_list[6]
        if strand == "+":
            orientation = "1"
        else:
            orientation = "-1"
        
        if re.search(r',similar to AA sequence:', line):
            
            db_name = re.search(',similar to AA sequence:(\w+)', line).group(1)
            
            if db_name == "user_db":
                gene = re.search(',similar to AA sequence:user_db:(\S+?);', line).group(1)
            elif db_name == "UniProtKB":
                if re.search(r';Name=', line):
                    gene = re.search(';Name=(\S+?);', line).group(1)
                    gene = gene.split('_')[0]
                else:
                    gene = re.search(',similar to AA sequence:UniProtKB:(\S+?);', line).group(1)
            
            else:
                gene = re.search(',similar to AA sequence:\w+?:(\S+?);', line).group(1)
                
            
            #test if overlap with abricate gene
            if (start < abricate_dict[molecule][0] and end < abricate_dict[molecule][0]) or (start > abricate_dict[molecule][1] and end > abricate_dict[molecule][1]):        
                outline = str(molecule) + "\t" + str(gene) + "\t" + str(start) + "\t" + str(end) + "\t" + str(orientation) + "\t" + str(db_name)
                print(outline, file=gggenes_file)
            
            
        elif re.search(r';product=hypothetical protein', line):
            
            gene = "Hypothetical protein"
            db_name = "None"
            
            #test if overlap with abricate gene
            if (start < abricate_dict[molecule][0] and end < abricate_dict[molecule][0]) or (start > abricate_dict[molecule][1] and end > abricate_dict[molecule][1]):
                outline = str(molecule) + "\t" + str(gene) + "\t" + str(start) + "\t" + str(end) + "\t" + str(orientation) + "\t" + str(db_name)
                print(outline, file=gggenes_file)

gff_file.close()
gggenes_file.close()
abricate_file.close()







