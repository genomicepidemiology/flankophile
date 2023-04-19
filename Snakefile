# FLANKOPHILE
# Alix Vincent Thorn

configfile: "config.yaml"

import os
import re

## Input control reference database ##############################################################

header_list = []
with open(config["database"], 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith(">"):
            if len(line) < 2:
                raise Exception("ERROR! Database headers must not be empty.")
            header_start = line[1]
            if not header_start.isalpha():
                if not header_start.isdigit(): 
                    raise Exception("ERROR! Database headers must start with a letter or a number.")
            if line.find(";") != -1:
                raise Exception("ERROR! Database headers must not contain semicolon.")
            if line.find(" ") != -1: # If whitespace present then everything after whitespace is header.
                line_list=line.split()
                this_header = line_list[0]
                if len(line_list[0]) > 60:
                    raise Exception("ERROR! Database headers must not be longer than 60 characters including >.")
            else: # if no whitespace in header
                this_header = line
                if len(line) > 60:
                    raise Exception("ERROR! Database headers must not be longer than 60 characters including >.")
            if this_header in header_list:
                error_uniq = "ERROR! Database headers must be unique. This header is not unique: " + str(this_header)
                raise Exception(error_uniq)
            else:
               header_list.append(this_header)
                
                


## Input control config file ################################################################

try:
  int(config["flank_length_upstreams"])
  int(config["flank_length_downstreams"])
except:
  print("ERROR! Flank length must be a positive integer.")
  


if int(config["flank_length_upstreams"]) == 0 and int(config["flank_length_downstreams"]) == 0:
    raise Exception("ERROR! Both flank_length_upstreams and flank_length_downstreams cannot both be 0. If you just want to find the genes and not the flanks use the tool Abricate by Torsten Seemann")

if int(config["flank_length_upstreams"]) < 0 or int(config["flank_length_downstreams"]) < 0:
    raise Exception("ERROR! Flank length must be a positive integer.")

try:
  float(config["min_coverage_abricate"])
except:
  print("ERROR! min_coverage_abricate must be a number, more specifically a float.")
  raise Exception("ERROR!")

if float(config["min_coverage_abricate"]) < 90 or float(config["min_coverage_abricate"]) > 100:
    raise Exception("ERROR! min_coverage_abricate must be between 90 and 100.")


try:
  float(config["min_identity_abricate"])
except:
  print("ERROR! min_identity_abricate must be a number, more specifically a float.")
  raise Exception("ERROR!")
    
if float(config["min_identity_abricate"]) < 90 or float(config["min_identity_abricate"]) > 100:
    raise Exception("ERROR! min_identity_abricate must be between 90 and 100.")


try:
  float(config["cluster_identity_and_length_diff"])  
except:
  print("ERROR! Cluster parameters must be numbers.")
  raise Exception("ERROR! cluster_identity_and_length_diff parameter must be a number." )

if float(config["cluster_identity_and_length_diff"]) < 0.8 or float(config["cluster_identity_and_length_diff"]) > 1:
    raise Exception("ERROR! cluster_identity_and_length_diff must be between 0.80 and 1.")
    
    

try:
  int(config["k-mer_size"])
  int(config["distance_measure"])
except:
  raise Exception("ERROR! k-mer_size and distance_measure must be a positive integer." )

if int(config["k-mer_size"]) < 6:
    raise Exception("ERROR! k-mer_size is wrong. K must be a positive integer above 5. Google k-mer size for guidance")

if int(config["distance_measure"]) == 1 or int(config["distance_measure"]) == 64 or int(config["distance_measure"]) == 256  or int(config["distance_measure"]) == 4096:
    status = "Good"
else:
    raise Exception("ERROR! distance_measure must be 1, 64, 256 or 4096. See README file.")

## Read in data #######################################################################################



# Decide on word size for kma

cd_id = float(config["cluster_identity_and_length_diff"])

if cd_id >= 0.95:
    word_size_cd_hit = 10
elif cd_id >= 0.90:
    word_size_cd_hit = 8
elif cd_id >= 0.88:
    word_size_cd_hit = 7
elif cd_id >= 0.85:
    word_size_cd_hit = 6
elif cd_id >= 0.80:
    word_size_cd_hit = 5
else:
    word_size_cd_hit = "None"
    



ASSEMBLY_NAMES = []
ASSEMBLY_NAME_PATH_DICT = {}
PATH_ASSEMBLY_NAME_DICT = {}

ASSEMBLY_NAME_META_DICT = {}

first_line_flag = False

col_num = "None"

with open(config["input_list"], 'r') as file:
    for line in file:
        line = line.strip()
        if not line.startswith("#"):
            line_list = line.split("\t")
            if first_line_flag == True:
                if len(line_list) != col_num:
                    raise Exception("ERROR! There must be the same number of columns on each row of input_list.")
            if len(line_list) > 3:
                raise Exception("ERROR! Input list must have 2 or 3 collums in each row and be tab separated. Too many collumns in this input file.")
            if len(line_list) < 2:
                raise Exception("ERROR! Input list must have 2 or 3 collums in each row and be tab separated. Too few collumns in this input file.")
            assembly_name = line_list[0]
            assembly_path = line_list[1]
            if assembly_path.find(" ") != -1:
                raise Exception("ERROR! Paths in the input_list file must not contain white space.")
            if assembly_name in ASSEMBLY_NAMES:
                raise Exception("ERROR! Assembly names have to be unique. All assembly names are not unique.")
            if assembly_name.find(" ") != -1:
                raise Exception("ERROR! Assembly names in the input_list must not contain white space.")
            if assembly_name.find("/") != -1:
                raise Exception("ERROR! Assembly names in the input_list must not contain slash.")
            
            ASSEMBLY_NAMES.append(assembly_name)
            ASSEMBLY_NAME_PATH_DICT[assembly_name] = assembly_path
            PATH_ASSEMBLY_NAME_DICT[assembly_path] = assembly_name
            first_line_flag = True
            col_num = len(line_list)
            if col_num == 3:
                metadata = line_list[2]
                meta_check = re.search('\W', metadata) # Check that metadata only contain letters, numbers and _
                if meta_check is None:
                    ASSEMBLY_NAME_META_DICT[assembly_name] = str(metadata)
                else: 
                    raise Exception("ERROR! Metadata column must only contain letters, numbers and underscore.")
            if col_num == 2:
                metadata = "No_metadata"
                
                
                

## Rule all ###########################################################################################

            
rule all:
    input: 
       "output/1_variants.fasta",
       "output/2_hits_included_in_flank_analysis.tsv",
       "output/3_clustering.tsv",
       "output/4_config.yaml"


## Flankophile script ##################################################################################

rule user_db:
    input:
        config["database"]
    output:
        temp("output/1_search/abricate.txt")
    conda: "environment.yaml"
    shell:
        "cat {input} | awk -F' ' '{{print $1}}' > bin/abricate/db/user_db/sequences;"
        "echo {input} > {output};"
        "./bin/abricate/bin/abricate --setupdb >> {output}"



rule abricate:
    input:
        db="output/1_search/abricate.txt"
    output:
        tsv=temp("output/1_search/tsv/{assembly_name}/{assembly_name}.tsv"),
        length="output/1_search/contig_length/{assembly_name}/{assembly_name}.length",
        tsv_c2=temp("output/1_search/tsv/{assembly_name}/{assembly_name}.tsv_c2"),
        length_temp=temp("output/1_search/contig_length/{assembly_name}/{assembly_name}.length_temp")
    params:
        assembly_path = lambda wildcards: ASSEMBLY_NAME_PATH_DICT[wildcards.assembly_name],
        cov=config["min_coverage_abricate"],
        id=config["min_identity_abricate"]
    conda: "environment.yaml"
    shell:
        "./bin/abricate/bin/abricate --db user_db --minid {params.id} --mincov {params.cov} {params.assembly_path} > {output.tsv};"
        "seqkit fx2tab --length --name  {params.assembly_path} | awk -v  OFS='\t' '{{print $1, $NF}}' > {output.length_temp};"
        "cut -f2 {output.tsv} > {output.tsv_c2};"
        "grep -wFf {output.tsv_c2} {output.length_temp} > {output.length} || true"





rule abricate_merge:
    input:
        expand("output/1_search/tsv/{assembly_name}/{assembly_name}.tsv", assembly_name=ASSEMBLY_NAMES)
    output:
         no_head=temp("output/1_search/abricate_all_raw_no_header.tsv"),
         tsv=temp("output/1_search/search_raw.tsv")
    shell:
        "touch {output.no_head};"
        "cat  {input} | cut -f1-11 | grep -v '^#'  >> {output.no_head} || true;"   # True to avoid error if no hits are found
        "cat bin/abricate_header.txt {output.no_head} > {output.tsv};"



rule python_enrich_abricate_output_with_length_and_metadata:
    input:
        "output/1_search/search_raw.tsv",
        expand("output/1_search/contig_length/{assembly_name}/{assembly_name}.length", assembly_name=ASSEMBLY_NAMES)
    output:
        temp("output/1_search/search_temp.tsv")            
    run:
        flag = False                                              # Flag remembers if any hits have been found
        input_tsv = open("output/1_search/search_raw.tsv", "r")
        output_tsv = open("output/1_search/search_temp.tsv","w") 
        for line in input_tsv:
            line = line = line.strip()
            if line.startswith("#"):
                print(line, file = output_tsv)
            else:
                line_list = line.split("\t")
                assembly_name = PATH_ASSEMBLY_NAME_DICT[line_list[0]]
                path_to_length_file = "output/1_search/contig_length/" + assembly_name + "/" + assembly_name + ".length"
                this_contig_len = "unknown"
                        
                input_length = open(path_to_length_file, "r")
                for line_l in input_length:
                    if this_contig_len == "unknown":
                        line_l = line_l.strip()
                        line_l_list = line_l.split("\t")           # Split the length file into contig name and length
                        if line_l_list[0] == line_list[1]:         # If it is the same contig name.
                            this_contig_len = int(line_l_list[1])  # Save the contig length.
                input_length.close
                
                if len(ASSEMBLY_NAME_META_DICT) == 0:
                    new_line = line + "\t" + str(this_contig_len) + "\t" + "No_metadata"
                else: 
                    new_line = line + "\t" + str(this_contig_len) + "\t" + ASSEMBLY_NAME_META_DICT[assembly_name]  
                print(new_line, file = output_tsv)
                flag = True                                        # Some hits were found.
        input_tsv.close
        output_tsv.close
        if flag == False:
            raise Exception("ERROR! ERROR! ERROR!\nERROR. Sorry. No hits found with requested min_coverage and min_identity values.")         
         


rule python_enrich_abricate_output_with_assembly_name_and_observation_ID:
    input:
        tsv="output/1_search/search_temp.tsv"
    output:
        tsv=temp("output/1_hits_all_no_var.tsv")
    run:
        input_tsv = open("output/1_search/search_temp.tsv", "r")
        output_tsv = open("output/1_hits_all_no_var.tsv","w")
        OBS_DICT = {}
        for line in input_tsv:
            line = line = line.strip()
            if line.startswith("#"):
                print(line, file = output_tsv)
            else:
                line_list = line.split("\t")
                assembly_name = PATH_ASSEMBLY_NAME_DICT[line_list[0]]
                
                if assembly_name[0:2].isalpha():
                    letter = assembly_name[0:2]
                elif assembly_name[0].isalpha():
                    letter = assembly_name[0]
                else:
                    letter = "i"
                    
                if letter in OBS_DICT:
                    OBS_DICT[letter] += 1
                else:
                    OBS_DICT[letter] = 1	

                new_line = line + "\t" + assembly_name + "\t" + letter + str(OBS_DICT[letter])
                print(new_line, file = output_tsv)
        input_tsv.close
        output_tsv.close
        






rule add_file_name_to_1_tsv:
    input:
        "output/1_hits_all_no_var.tsv"
    output:
        temp("output/temp_variant/all_hits_with_future_filename.tsv")
    run:
        input=open(input[0], "r")
        output=open(output[0], "w")
        for line in input:
            line = line.strip()
            if not line.startswith("#"):
                flag = True
                line_list = line.split("\t")
                ASSEMBLY_NAME = line_list[13]
                future_file_name = ASSEMBLY_NAME + ".tsv"
                new_line = line + "\t" + future_file_name
                print(new_line, file=output) 
        input.close()
        output.close()


checkpoint split_abricate_results1:
    input:
        "output/temp_variant/all_hits_with_future_filename.tsv"
    output:
        clusters=directory("output/temp_variant/abricate_results_per_assembly")
    shell:
        "mkdir output/temp_variant/abricate_results_per_assembly;"
        "cd output/temp_variant/abricate_results_per_assembly;"
        "cat ../all_hits_with_future_filename.tsv  | awk  '{{print>$16}}'"
        
        



rule make_bedfiles1:
    input:
        "output/temp_variant/abricate_results_per_assembly/{a}.tsv"
    output:
        temp("output/temp_variant/bedfiles/target_sequence_only/{a}/{a}.bed")
    run:
        input=open(input[0], "r")
        output_target_sequence_only=open(output[0], "w")
        col_index_SEQ, col_index_START, col_index_END, col_index_STRAND, col_index_ID_num  = 1, 2, 3, 4, 14

        for line in input:
            line = line.strip()
            line_list = line.split()
            contig_name = line_list[col_index_SEQ]
            START_old = int(line_list[col_index_START])
            END_old = int(line_list[col_index_END])
            gene_strand = line_list[col_index_STRAND]
            ID = line_list[col_index_ID_num]
            # cut out the gene  with bedtools getfasta
            out_line = contig_name +  "\t" + str(START_old - 1) +  "\t" + str(END_old)
            out_line = out_line + "\t" + ID + "\t" + "1" + "\t" + gene_strand
            print(out_line, file=output_target_sequence_only)
        input.close()
        output_target_sequence_only.close





rule bedtools_target_sequence_only1:
    input:
        bed="output/temp_variant/bedfiles/target_sequence_only/{a}/{a}.bed"
    output:
        fasta=temp("output/temp_variant/gene_fasta_files_per_assembly/target_sequence_only/{a}/{a}.fa")
    conda: "environment.yaml"
    params:
        assembly_path=lambda wildcards: ASSEMBLY_NAME_PATH_DICT[wildcards.a]
    shell:
        '''
        bedtools getfasta -fi {params.assembly_path} -bed {input.bed} -s -nameOnly > {output.fasta};
        sed -i -e "s/[(|)|+]//g" {output.fasta} && sed -i -e "s/-//g" {output.fasta}
        ''' 



def target_sequence_only_input1(wildcards):
    checkpoint_output = checkpoints.split_abricate_results1.get(**wildcards).output[0]
    return expand("output/temp_variant/gene_fasta_files_per_assembly/target_sequence_only/{a}/{a}.fa",
           a=glob_wildcards(os.path.join(checkpoint_output, "{a}.tsv")).a)


rule pool_target_sequence_only1:
    input:
        target_sequence_only_input1
    output:
        temp=temp("output/temp_variant/gene_fasta_files_pooled/temp_target_sequence_only.fa"),
        fasta=temp("output/temp_variant/target_sequence_only.fa")
    conda: "environment.yaml"
    shell:
        "rm -rf output/temp_variant/abricate_results_per_assembly;" 
        "cat {input}  > {output.temp};"
        "cat {output.temp} | seqkit sort -N -i -2  > {output.fasta}"


rule merge_duplicate_sequences1:
    input:
        "output/temp_variant/target_sequence_only.fa"
    output:
        temp("output/temp_variant/target_sequence_only_no_dup.fa")
    shell:
        "python3 bin/DupRemover.py -i {input} -o {output};"
        
        
rule add_variant_numbers_to_fasta1:
    input:
        "output/temp_variant/target_sequence_only_no_dup.fa"
    output:
        temp("output/temp_variant/target_sequence_only_no_dup_var_num.fa")
    run:
        input_file=open(input[0], "r")
        output_file=open(output[0], "w")
        counter_fasta = 0
        for line in input_file:
            line = line.strip()
            if line.startswith(">"):
                counter_fasta += 1
                novel_line = line + "|v_" + str(counter_fasta)
                print(novel_line, file = output_file)
            else:
                print(line, file = output_file)
        input_file.close()
        output_file.close()





rule make_variant_fasta:
    input:
        "output/temp_variant/target_sequence_only_no_dup_var_num.fa"
    output:
        "output/1_variants.fasta"
    run:
        input_file=open(input[0], "r")
        output_file=open(output[0], "w")
        for line in input_file:
            line = line.strip()
            if line.startswith(">"):
                line = line[1:]                    # Remove >
                line_list = line.split("|")
                name_of_variant = line_list[-1]
                new_header_line = ">" + name_of_variant
                print(new_header_line, file = output_file)
            else:
                print(line, file = output_file)
        input_file.close()
        input_file.close()
        
        
        
        

rule add_variant_number_to_results1:
    input:
        fasta="output/temp_variant/target_sequence_only_no_dup_var_num.fa",
        tsv="output/1_hits_all_no_var.tsv"
    output:
        tsv="output/1_hits_all.tsv"
    run:
        input_tsv=open("output/1_hits_all_no_var.tsv", "r")
        output_tsv=open("output/1_hits_all.tsv", "w")
        
        for tsv_line in input_tsv:
            variants_name = "None"
            tsv_line = tsv_line.strip()
            if tsv_line.startswith("#"):
                print(tsv_line, file = output_tsv)
            else:
                tsv_line_list = tsv_line.split("\t")
                obs_id = tsv_line_list[14]
                input_fasta=open("output/temp_variant/target_sequence_only_no_dup_var_num.fa", "r")
                for fasta_line in input_fasta:                   
                    if fasta_line.startswith(">"):
                        fasta_line = fasta_line.strip()[1:]     #remove the >
                        fasta_line_list = fasta_line.split("|")
                        if obs_id in fasta_line_list:           # if it is the correct variant
                            variants_name = fasta_line_list[-1]
                new_tsv_line = tsv_line + "\t" + variants_name 
                print(new_tsv_line, file = output_tsv)
                input_fasta.close()
        input_tsv.close()
        output_tsv.close()
        
        
        

rule filter_abricate_space_for_flanks:
    input:
        tsv="output/1_hits_all.tsv"
    output:
        tsv="output/2_hits_included_in_flank_analysis.tsv",
        report="output/2_report_flank_filtering.txt"
    run:
        counter_accepted = 0

        # Counter of how many observations that would be discarted at each threshold
        UP_user, UP_half, UP_double, DOWN_user, DOWN_half, DOWN_double = 0, 0, 0, 0, 0, 0 
        
        total_input_observations, total_discarded_observation = 0, 0  
        col_index_ASSEMBLY_PATH, col_index_SEQ, col_index_START, col_index_END, col_index_STRAND  = 0, 1, 2, 3, 4
        user_upstreams=float(config["flank_length_upstreams"])
        user_downstreams=float(config["flank_length_downstreams"])
                 

        input_tsv = open("output/1_hits_all.tsv", "r")
        output_tsv = open("output/2_hits_included_in_flank_analysis.tsv", "w")
        for line in input_tsv:
            line = line.strip()
            if line.startswith("#"): print(line, file = output_tsv)
            else:
                space_for_flank_up, space_for_flank_down = 0, 0
                total_input_observations += 1
                line_list = line.split("\t")
                path = line_list[col_index_ASSEMBLY_PATH]
                seq_name = line_list[col_index_SEQ]
                gene_start = int(line_list[col_index_START])
                gene_end = int(line_list[col_index_END])
                gene_strand = line_list[col_index_STRAND]
                this_contig_len = int(line_list[11])                #CONTIG_LENGTH column in tsv
                                                                                  
                if gene_strand == "+":
                    space_for_flank_up = gene_start - 1
                    space_for_flank_down = this_contig_len - gene_end
                if gene_strand == "-":
                    space_for_flank_down = gene_start - 1
                    space_for_flank_up = this_contig_len - gene_end
                if space_for_flank_up < round(user_upstreams * 2):
                    UP_double += 1
                    if space_for_flank_up < user_upstreams:
                        UP_user += 1
                        if space_for_flank_up < round(user_upstreams / 2):
                            UP_half += 1
                if space_for_flank_down < round(user_downstreams * 2):
                    DOWN_double += 1
                    if space_for_flank_down < user_downstreams:
                        DOWN_user += 1
                        if space_for_flank_down < round(user_downstreams / 2):
                            DOWN_half += 1
                if  space_for_flank_up >= user_upstreams and space_for_flank_down >= user_downstreams:
                    counter_accepted += 1

                    print(line, file = output_tsv)
                else:
                    total_discarded_observation += 1
        input_tsv.close
        output_tsv.close
        output_report = open("output/2_report_flank_filtering.txt", "w")
        a = "Flank length upstreams:" + "\t" + str(int(user_upstreams))
        b = a +"\nFlank length downstreams:" + "\t" + str(int(user_downstreams))
        c = b + "\n\nNumber of hits before flank space filtering:" + "\t" + str(total_input_observations)
        d = c + "\nNumber of hits after flank space filtering:" + "\t" + str(total_input_observations - total_discarded_observation)
        e = d + "\n\nNumber of hits discarded by flank space filtering:" + "\t" + str(total_discarded_observation)
        f = e +"\n\n\nHypothetical number of hits discarded due to insufficient upstreams flank space if the following thresholds were used:"
        g = f + "\n" + str(round(user_upstreams / 2)) + "\t" + str(round(user_upstreams)) + "\t"  + str(round(user_upstreams * 2)) 
        h = g + "\n" +  str(UP_half) + "\t" + str(UP_user) + "\t"  + str(UP_double)
        i = h + "\n\nHypothetical number of hits discarded due to insufficient downstreams flank space if the following thresholds were used:"
        j = i + "\n" + str(round(user_downstreams / 2)) + "\t" + str(round(user_downstreams)) + "\t"  + str(round(user_downstreams* 2))
        k = j + "\n" +  str(DOWN_half) + "\t" + str(DOWN_user) + "\t"  + str(DOWN_double)
        print(k, file = output_report)
        output_report.close()



rule add_file_name_to_2_tsv:
    input:
        "output/2_hits_included_in_flank_analysis.tsv"
    output:
        temp("output/temp_1/final_gene_results_with_future_filename.tsv")
    run:
        flag = False
        input=open(input[0], "r")
        output=open(output[0], "w")
        for line in input:
            line = line.strip()
            if not line.startswith("#"):
                flag = True
                line_list = line.split()
                ASSEMBLY_NAME = line_list[13]
                future_file_name = ASSEMBLY_NAME + ".tsv"
                new_line = line + "\t" + future_file_name
                print(new_line, file=output) 
        input.close()
        output.close()
        if flag == False:
            raise Exception("ERROR! ERROR! ERROR!\nERROR. Sorry. No hits found with requested flanking regions length. Delete output folder and try again with a lower flank length.")
            
        
        
checkpoint split_abricate_results2:
    input:
        "output/temp_1/final_gene_results_with_future_filename.tsv"
    output:
        clusters=directory("output/temp_1/1_abricate_results_per_assembly")
    shell:
        "rm -rf output/1_search;"
        "mkdir output/temp_1/1_abricate_results_per_assembly;"
        "cd output/temp_1/1_abricate_results_per_assembly;"
        "cat ../final_gene_results_with_future_filename.tsv  | awk  '{{print>$17}}'"
        


rule make_bedfiles2:
    input:
        "output/temp_1/1_abricate_results_per_assembly/{a}.tsv"
    output:
        temp("output/temp_2/bedfiles/target_and_flanking_regions/{a}/{a}.bed"),
        temp("output/temp_2/bedfiles/target_sequence_only/{a}/{a}.bed"),
        temp("output/temp_2/bedfiles/flanking_regions_only/{a}/{a}.bed")
    run:
        input=open(input[0], "r")
        output_target_and_flanking_regions=open(output[0], "w")
        output_target_sequence_only=open(output[1], "w")
        output_flanking_regions_only=open(output[2], "w")
        col_index_SEQ, col_index_START, col_index_END, col_index_STRAND, col_index_ID_num  = 1, 2, 3, 4, 14
        user_upstreams=int(config["flank_length_upstreams"])
        user_downstreams=int(config["flank_length_downstreams"])

        for line in input:
            line = line.strip()
            line_list = line.split()
            contig_name = line_list[col_index_SEQ]
            START_old = int(line_list[col_index_START])
            END_old = int(line_list[col_index_END])
            gene_strand = line_list[col_index_STRAND]
            ID = line_list[col_index_ID_num]

            # for target_and_flanking_regions - whole region - with bedtools getfasta
            if gene_strand == "+":
                START_target_and_flanking_regions = START_old - user_upstreams -1
                END_target_and_flanking_regions = END_old + user_downstreams
            if gene_strand == "-":
                START_target_and_flanking_regions = START_old - user_downstreams -1
                END_target_and_flanking_regions = END_old + user_upstreams
            out_line = contig_name +  "\t" + str(START_target_and_flanking_regions) +  "\t" + str(END_target_and_flanking_regions)
            out_line = out_line + "\t" + ID + "\t" + "1" + "\t" + gene_strand
            print(out_line, file=output_target_and_flanking_regions)
 
            #target_sequence_only - just cut out the gene  with bedtools getfasta
            out_line = contig_name +  "\t" + str(START_old - 1) +  "\t" + str(END_old)
            out_line = out_line + "\t" + ID + "\t" + "1" + "\t" + gene_strand
            print(out_line, file=output_target_sequence_only)

            #flanking_regions_only - made with maskfasta from target_and_flanking_regions_fasta. All strands will be positive
            START_flanking_regions_only = user_upstreams + 1 - 1
            gene_length = END_old - START_old + 1
            END_flanking_regions_only = user_upstreams + gene_length
            out_line = ID +  "\t" + str(START_flanking_regions_only) + "\t" + str(END_flanking_regions_only)
            print(out_line, file =output_flanking_regions_only )
        input.close()
        output_target_and_flanking_regions.close
        output_target_sequence_only.close
        output_flanking_regions_only.close
        


rule bedtools_target_and_flanking_regions:
    input:
        bed="output/temp_2/bedfiles/target_and_flanking_regions/{a}/{a}.bed"
    output:
        fasta=temp("output/temp_2/gene_fasta_files_per_assembly/target_and_flanking_regions/{a}/{a}.fa")
    conda: "environment.yaml"
    params:
        assembly_path=lambda wildcards: ASSEMBLY_NAME_PATH_DICT[wildcards.a]
    shell:
        '''
        bedtools getfasta -fi {params.assembly_path} -bed {input.bed} -s -nameOnly > {output.fasta};
        sed -i -e "s/[(|)|+]//g" {output.fasta} && sed -i -e "s/-//g" {output.fasta}
        '''

rule bedtools_target_sequence_only2:
    input:
        whole_flank="output/2_filter/gene_fasta_files_pooled/target_and_flanking_regions.fa",
        bed="output/temp_2/bedfiles/target_sequence_only/{a}/{a}.bed"
    output:
        fasta=temp("output/temp_2/gene_fasta_files_per_assembly/target_sequence_only/{a}/{a}.fa")
    conda: "environment.yaml"
    params:
        assembly_path=lambda wildcards: ASSEMBLY_NAME_PATH_DICT[wildcards.a]
    shell:
        '''
        bedtools getfasta -fi {params.assembly_path} -bed {input.bed} -s -nameOnly > {output.fasta};
        sed -i -e "s/[(|)|+]//g" {output.fasta} && sed -i -e "s/-//g" {output.fasta}
        ''' 


rule bedtools_flanking_regions_only:
    input:
        target_sequence_only="output/2_filter/gene_fasta_files_pooled/target_sequence_only.fa",
        bed="output/temp_2/bedfiles/flanking_regions_only/{a}/{a}.bed",
        fasta="output/temp_2/gene_fasta_files_per_assembly/target_and_flanking_regions/{a}/{a}.fa"
    output:
        fasta=temp("output/temp_2/gene_fasta_files_per_assembly/flanking_regions_only/{a}/{a}.fa")
    conda: "environment.yaml"
    shell:
        "bedtools maskfasta -fi {input.fasta} -bed {input.bed} -fo {output.fasta}"



def target_and_flanking_regions_input(wildcards):
    checkpoint_output = checkpoints.split_abricate_results2.get(**wildcards).output[0]
    return expand("output/temp_2/gene_fasta_files_per_assembly/target_and_flanking_regions/{a}/{a}.fa",
           a=glob_wildcards(os.path.join(checkpoint_output, "{a}.tsv")).a)



def target_sequence_only_input(wildcards):
    checkpoint_output = checkpoints.split_abricate_results2.get(**wildcards).output[0]
    return expand("output/temp_2/gene_fasta_files_per_assembly/target_sequence_only/{a}/{a}.fa",
           a=glob_wildcards(os.path.join(checkpoint_output, "{a}.tsv")).a)


def flanking_regions_only_input(wildcards):
    checkpoint_output = checkpoints.split_abricate_results2.get(**wildcards).output[0]
    return expand("output/temp_2/gene_fasta_files_per_assembly/flanking_regions_only/{a}/{a}.fa",
           a=glob_wildcards(os.path.join(checkpoint_output, "{a}.tsv")).a)


rule pool_target_and_flanking_regions:
    input:
        target_and_flanking_regions_input
    output:
        temp=temp("output/temp_2/gene_fasta_files_pooled/temp_target_and_flanking_regions.fa"),
        fasta="output/2_filter/gene_fasta_files_pooled/target_and_flanking_regions.fa"
    conda: "environment.yaml"
    shell:
        "cat {input}  > {output.temp};"
        "cat {output.temp} | seqkit sort -N -i -2  > {output.fasta};"
        "rm -r output/temp_1"


rule pool_target_sequence_only2:
    input:
        target_sequence_only_input
    output:
        temp=temp("output/temp_2/gene_fasta_files_pooled/temp_target_sequence_only.fa"),
        fasta="output/2_filter/gene_fasta_files_pooled/target_sequence_only.fa"
    conda: "environment.yaml"
    shell:
        "cat {input}  > {output.temp};"
        "cat {output.temp} | seqkit sort -N -i -2  > {output.fasta}"


rule pool_flanking_regions_only:
    input:
        flanking_regions_only_input
    output:
        temp=temp("output/temp_2/gene_fasta_files_pooled/temp_flanking_regions_only.fa"),
        fasta="output/2_filter/gene_fasta_files_pooled/flanking_regions_only.fa"
    conda: "environment.yaml"
    shell:
        "cat {input}  > {output.temp};"
        "cat {output.temp} | seqkit sort -N -i -2  > {output.fasta}"


rule extract_relavant_reference_sequences:
    input:
        flanks="output/2_filter/gene_fasta_files_pooled/target_and_flanking_regions.fa",
        target_sequence_only="output/2_filter/gene_fasta_files_pooled/target_sequence_only.fa",
        masked="output/2_filter/gene_fasta_files_pooled/flanking_regions_only.fa",
        final_results="output/2_hits_included_in_flank_analysis.tsv"
    output:
        list=temp("output/3_define_clusters/final_output_list.txt"),
        fasta=temp("output/3_define_clusters/relevant_ref_seqs.fasta")
    conda: "environment.yaml"
    shell:
        "awk '{{print $6}}' {input.final_results} |awk 'NR>1' | sort -u > {output.list};"
        "seqkit grep -n -f {output.list} bin/abricate/db/user_db/sequences > {output.fasta}"


rule cd_hit:
    input:
        "output/3_define_clusters/relevant_ref_seqs.fasta"
    output:
        clus=temp("output/3_define_clusters/cd_hit.clstr"),
        genes=temp("output/3_define_clusters/cd_hit"),
        new_head=temp("output/3_define_clusters/cd_hit.txt"),
        temp=temp("output/3_define_clusters/cd_hit.temp")
    conda: "environment.yaml"
    params:
        c_s= config["cluster_identity_and_length_diff"],
        n = word_size_cd_hit
    shell:
        '''
        cd-hit-est -i {input} -o {output.genes} -M 0 -d 0 -c {params.c_s} -n {params.n} -s {params.c_s}  -sc 1 -g 1;
        perl bin/clstr2txt.pl {output.clus} | tail -n+2 > {output.temp};
        cat bin/cd-hit-txt_header.txt {output.temp} > {output.new_head}
        '''


rule process_cd_hit_results:
    input: 
        "output/3_define_clusters/cd_hit.txt"
    output:
        "output/3_clustering.tsv"
    run:
        dict = {}
        input=open(input[0], "r")
        output=open(output[0], "w")

        for line in input:
            if line.startswith("#"):
                line = line.strip()
                print(line, file=output)
            if not line.startswith("#"):
                line = line.strip()
                line_list = line.split()
                gene_name = line_list[0]
                clus_number = line_list[1]
                if clus_number not in dict:
                    # Transform characters that are not file name friendly into _                    
                    friendly_gene_name = re.sub(r'\W', '_', gene_name)
                    dict[clus_number] = clus_number + "_" + friendly_gene_name 
                new_line = line + "\t" + dict[clus_number] + ".tsv" + "\t" + dict[clus_number]
                print(new_line, file=output)        
        input.close()
        output.close()


checkpoint split_cd_hit_results:
    input:
        "output/3_clustering.tsv"
    output:
        clusters=directory("output/3_define_clusters/cd_hit_per_cluster")
    shell:
        "mkdir output/3_define_clusters/cd_hit_per_cluster;"
        "cd output/3_define_clusters/cd_hit_per_cluster;"
        "cat ../../3_clustering.tsv | grep -v '^#'| awk  '{{print>$8}}'"




rule extract_cluster_rows_from_results_file:
    input:
        results="output/2_hits_included_in_flank_analysis.tsv",
        list="output/3_define_clusters/cd_hit_per_cluster/{c}.tsv"
    output:
        results="output/4_cluster_results/{c}/{c}.tsv",
        temp=temp("output/4_cluster_results/{c}/{c}.temp")
    params:
        awk="'NR==FNR{a[$1];next}($6 in a){{print $0}}'"        # 6 refer to GENE collumn in "output/2_hits_included_in_flank_analysis.tsv"
    shell:
        '''
        awk {params.awk} {input.list} {input.results} > {output.temp};
        cat bin/abricate_header.txt {output.temp} > {output.results}
        '''


rule extract_cluster_fastas:
    input:
        tsv="output/4_cluster_results/{c}/{c}.tsv",
        target_and_flanking_regions_fasta="output/2_filter/gene_fasta_files_pooled/target_and_flanking_regions.fa",
        target_sequence_only_fasta="output/2_filter/gene_fasta_files_pooled/target_sequence_only.fa",
        flanking_regions_only_fasta="output/2_filter/gene_fasta_files_pooled/flanking_regions_only.fa"
    output:
        ID_list=temp("output/4_cluster_results/{c}/{c}.ID_list"),
        target_and_flanking_regions_fasta="output/4_cluster_results/{c}/{c}.target_and_flanking_regions_fasta",
        target_sequence_only_fasta="output/4_cluster_results/{c}/{c}.target_sequence_only_fasta",
        flanking_regions_only_fasta="output/4_cluster_results/{c}/{c}.flanking_regions_only_fasta"
    conda: "environment.yaml"
    shell:
        "awk '{{print $15}}' {input.tsv} > {output.ID_list};"
        "seqkit grep -n -f {output.ID_list} {input.target_and_flanking_regions_fasta} > {output.target_and_flanking_regions_fasta};"
        "seqkit grep -n -f {output.ID_list} {input.target_sequence_only_fasta} > {output.target_sequence_only_fasta};"
        "seqkit grep -n -f {output.ID_list} {input.flanking_regions_only_fasta} > {output.flanking_regions_only_fasta}"


  

def split_target_sequence_only_fasta_file_input(wildcards):
    checkpoint_output = checkpoints.split_cd_hit_results.get(**wildcards).output[0]
    return expand("output/4_cluster_results/{c}/{c}.target_sequence_only_fasta",
           c=glob_wildcards(os.path.join(checkpoint_output, "{c}.tsv")).c)


rule translate_reference_database:
    input:
        split_target_sequence_only_fasta_file_input
    output:
        temp("output/user_db")
    conda: "environment.yaml"
    shell:
        "seqkit translate bin/abricate/db/user_db/sequences > {output}"


rule prokka:
    input:
        AA_db="output/user_db",
        fasta="output/4_cluster_results/{c}/{c}.target_and_flanking_regions_fasta"
    output:
        gff="output/4_cluster_results/{c}/prokka/{c}.gff",
        gbk="output/4_cluster_results/{c}/prokka/{c}.gbk",
        p1=temp("output/4_cluster_results/{c}/prokka/proteins.pot"),
        p2=temp("output/4_cluster_results/{c}/prokka/proteins.pdb"),
        p3=temp("output/4_cluster_results/{c}/prokka/proteins.ptf"),
        p4=temp("output/4_cluster_results/{c}/prokka/proteins.pto"),
        err=temp("output/4_cluster_results/{c}/prokka/{c}.err"),
        log=temp("output/4_cluster_results/{c}/prokka/{c}.log"),
        faa=temp("output/4_cluster_results/{c}/prokka/{c}.faa"),
        fnn=temp("output/4_cluster_results/{c}/prokka/{c}.ffn"),
        fna=temp("output/4_cluster_results/{c}/prokka/{c}.fna"),
        txt=temp("output/4_cluster_results/{c}/prokka/{c}.txt"),
        tbl=temp("output/4_cluster_results/{c}/prokka/{c}.tbl"),
        sqn=temp("output/4_cluster_results/{c}/prokka/{c}.sqn"),
        tsv=temp("output/4_cluster_results/{c}/prokka/{c}.tsv"),
        fsa=temp("output/4_cluster_results/{c}/prokka/{c}.fsa")
    params:
        outdir=lambda wildcards: "output/4_cluster_results/" +  wildcards.c + "/prokka",
        c=lambda wildcards: wildcards.c
    conda: "environment.yaml" 
    shell:
        "prokka --proteins {input.AA_db} --outdir {params.outdir} --force --prefix {params.c} {input.fasta}"



rule make_gggenes_input:
    input:
        gff="output/4_cluster_results/{c}/prokka/{c}.gff",
        tsv="output/4_cluster_results/{c}/{c}.tsv"
    output:
        g="output/4_cluster_results/{c}/{c}.gggenes"
    shell:
        "python3 ./bin/gff_to_gggenes.py {input.gff} {output.g} {input.tsv}"



rule kma_index:
    input:
        gggenes="output/4_cluster_results/{c}/{c}.gggenes",
        target_and_flanking_regions="output/4_cluster_results/{c}/{c}.target_and_flanking_regions_fasta",
        masked="output/4_cluster_results/{c}/{c}.flanking_regions_only_fasta",
        target_sequence_only="output/4_cluster_results/{c}/{c}.target_sequence_only_fasta"
    output:
        temp("output/4_cluster_results/{c}/kma_index/target_and_flanking_regions.comp.b"),
        temp("output/4_cluster_results/{c}/kma_index/target_and_flanking_regions.length.b"),
        temp("output/4_cluster_results/{c}/kma_index/target_and_flanking_regions.name"),
        temp("output/4_cluster_results/{c}/kma_index/target_and_flanking_regions.seq.b"),
        temp("output/4_cluster_results/{c}/kma_index/flanking_regions_only.comp.b"),
        temp("output/4_cluster_results/{c}/kma_index/flanking_regions_only.length.b"),
        temp("output/4_cluster_results/{c}/kma_index/flanking_regions_only.name"),
        temp("output/4_cluster_results/{c}/kma_index/flanking_regions_only.seq.b"),
        temp("output/4_cluster_results/{c}/kma_index/target_sequence_only.comp.b"),
        temp("output/4_cluster_results/{c}/kma_index/target_sequence_only.length.b"),
        temp("output/4_cluster_results/{c}/kma_index/target_sequence_only.name"),
        temp("output/4_cluster_results/{c}/kma_index/target_sequence_only.seq.b")
    conda: "environment.yaml"
    params:
        outfolder_target_and_flanking_regions=lambda wildcards: "output/4_cluster_results/" +  wildcards.c + "/kma_index/target_and_flanking_regions",
        outfolder_masked=lambda wildcards: "output/4_cluster_results/" +  wildcards.c + "/kma_index/flanking_regions_only",
        outfolder_target_sequence_only=lambda wildcards: "output/4_cluster_results/" +  wildcards.c + "/kma_index/target_sequence_only",
        k=config["k-mer_size"]
    shell:
        "kma index -k {params.k} -i {input.target_and_flanking_regions} -o {params.outfolder_target_and_flanking_regions};"
        "kma index -k {params.k} -i {input.masked} -o {params.outfolder_masked};"
        "kma index -k {params.k} -i {input.target_sequence_only} -o {params.outfolder_target_sequence_only}"




rule kma_dist:
    input:
        ja="output/4_cluster_results/{c}/kma_index/target_sequence_only.comp.b",
        jb="output/4_cluster_results/{c}/kma_index/target_sequence_only.length.b",
        jc="output/4_cluster_results/{c}/kma_index/target_sequence_only.name",
        fa="output/4_cluster_results/{c}/kma_index/target_and_flanking_regions.comp.b",
        fb="output/4_cluster_results/{c}/kma_index/target_and_flanking_regions.length.b",
        fc="output/4_cluster_results/{c}/kma_index/target_and_flanking_regions.name",
        ma="output/4_cluster_results/{c}/kma_index/flanking_regions_only.comp.b",
        mb="output/4_cluster_results/{c}/kma_index/flanking_regions_only.length.b",
        mc="output/4_cluster_results/{c}/kma_index/flanking_regions_only.name",
        f="output/4_cluster_results/{c}/kma_index/target_and_flanking_regions.seq.b",
        m="output/4_cluster_results/{c}/kma_index/flanking_regions_only.seq.b",
        j="output/4_cluster_results/{c}/kma_index/target_sequence_only.seq.b"
    output:
        f="output/4_cluster_results/{c}/{c}.target_and_flanking_regions_dist",
        m="output/4_cluster_results/{c}/{c}.flanking_regions_only_dist",
        j="output/4_cluster_results/{c}/{c}.target_sequence_only_dist"
    params:
        f_db=lambda wildcards: "output/4_cluster_results/" +  wildcards.c + "/kma_index/target_and_flanking_regions",
        m_db=lambda wildcards: "output/4_cluster_results/" +  wildcards.c + "/kma_index/flanking_regions_only",
        j_db=lambda wildcards: "output/4_cluster_results/" +  wildcards.c + "/kma_index/target_sequence_only",
        dist_measure=config["distance_measure"]
     
    conda: "environment.yaml"
    shell:
        "kma dist -t_db {params.f_db} -o {output.f} -d {params.dist_measure};"
        "kma dist -t_db {params.m_db} -o {output.m} -d {params.dist_measure};"
        "kma dist -t_db {params.j_db} -o {output.j} -d {params.dist_measure}"




def dist_input(wildcards):
    checkpoint_output = checkpoints.split_cd_hit_results.get(**wildcards).output[0]
    return expand("output/4_cluster_results/{c}/{c}.flanking_regions_only_dist",
           c=glob_wildcards(os.path.join(checkpoint_output, "{c}.tsv")).c)



rule plots:
    input:
        dist_input
    output:
        config="output/4_config.yaml",
        txt="output/4_plots/version.txt"
    conda: "environment.yaml"
    shell:
        '''
        cp config.yaml {output.config};
        rm -r output/3_define_clusters/cd_hit_per_cluster;
        rm -rf output/2_filter;
        cp {output.config} {output.txt};
        Rscript bin/plot_gene_clusters_from_flankophile.R
        '''


heart = "\n  #####   #####\n ####### #######\n#################\n## Flankophile ##\n ###############\n  #############"
heart = heart + "\n   ###########\n    #########\n     #######\n      #####\n       ###\n        #"

onsuccess:
    print("Flankophile finished. Thank you for using Flankophile.")
    print(heart)

