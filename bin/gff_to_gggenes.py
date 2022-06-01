# By Alex Vincent Thorn


import os

import sys

import argparse

import re

gff_path = str(sys.argv[1])

gggenes_path = str(sys.argv[2])



gff_file = open(gff_path,"r")

gggenes_file = open(gggenes_path,"w")

outheader = "molecule" + "\t" + "gene" + "\t" + "start" + "\t" + "end" + "\t" + "orientation" + "\t" + "db_name"

print(outheader, file=gggenes_file)

fasta_flag = False

for line in gff_file:
    line = line.strip()
    if line.startswith("##FASTA"):
        fasta_flag = True
    if fasta_flag == False and not line.startswith("##"):

        
        if re.search(r',similar to AA sequence:', line):
            line_list = line.split()
            molecule = line_list[0]
            start = line_list[3]
            end = line_list[4]
            strand = line_list[6]
            if strand == "+":
                orientation = "1"
            else:
                orientation = "-1"
            
            db_name = re.search(',similar to AA sequence:(\w+)', line).group(1)
            
            if db_name == "user_db":
                gene = re.search(',similar to AA sequence:user_db:(\S+?);', line).group(1)
            elif db_name == "UniProtKB":
                if re.search(r';Name=', line):
                    gene = re.search(';Name=(\w+?);', line).group(1)
                    gene = gene.split('_')[0]
                else:
                    gene = re.search(',similar to AA sequence:UniProtKB:(\S+?);', line).group(1)
            
            else:
                gene = re.search(',similar to AA sequence:\w+?:(\S+?);', line).group(1)
            
            outline = molecule + "\t" + gene + "\t" + start + "\t" + end + "\t" + orientation + "\t" + db_name
            print(outline, file=gggenes_file)
        

gff_file.close()
gggenes_file.close()



