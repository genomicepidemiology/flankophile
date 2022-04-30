# FLANKOPHILE version 0.0.2
# Alex Thorn

configfile: "config.yaml"

import os

ASSEMBLY_NAMES = []
ASSEMBLY_NAME_PATH_DICT = {}
PATH_ASSEMBLY_NAME_DICT = {}
ASSEMBLY_NAME_CONTIG_LIST_DICT = {}

with open(config["input_list"], 'r') as file:
    for line in file:
        line = line.strip()
        if not line.startswith("#"):
            line_list = line.split()
            assembly_name = line_list[0]
            assembly_path = line_list[-1]
            if assembly_name not in ASSEMBLY_NAMES:
                ASSEMBLY_NAMES.append(assembly_name)
                ASSEMBLY_NAME_PATH_DICT[assembly_name] = assembly_path
                PATH_ASSEMBLY_NAME_DICT[assembly_path] = assembly_name
            if config["input_format"] == "contigs":
                contig_name = line_list[1]
                if assembly_name in ASSEMBLY_NAME_CONTIG_LIST_DICT:
                    ASSEMBLY_NAME_CONTIG_LIST_DICT[assembly_name].append(contig_name)
                else: ASSEMBLY_NAME_CONTIG_LIST_DICT[assembly_name] = [contig_name]
            
rule all:
    input: 
       # Run abricate for all samples and filter results
       "output/2_filter_gene_observations/3_final_gene_results.tsv",
       # extract flanks and gene sequences and cluster refernce genes by identity.
       "output/4_cluster_by_gene_family/cd_hit.processed",
       # Make gene annotation and distance trees
       "output/99_trees_and_distance_matrixes_done",
       # Last step, always run. Save config file in output folder so record is kept
       "output/99_config.yaml"


####################################################################################

rule user_db:
    input:
        config["database"]
    output:
        "output/0_setup_abricate_db/abricate.txt"
    conda: "environment.yaml"
    shell:
        "cp {input} bin/abricate/db/user_db/sequences;"
        "echo {input} > {output};"
        "./bin/abricate/bin/abricate --setupdb >> {output}"




rule abricate_a:
    input:
        db="output/0_setup_abricate_db/abricate.txt"
    output:
        tsv="output/1_abricate_a/{assembly_name}/{assembly_name}.tsv",
        length="output/1_abricate_a/{assembly_name}/{assembly_name}.length"
    params:
        assembly_path = lambda wildcards: ASSEMBLY_NAME_PATH_DICT[wildcards.assembly_name]
    conda: "environment.yaml"
    shell:
        "./bin/abricate/bin/abricate --db user_db {params.assembly_path} > {output.tsv};"
        "seqkit fx2tab --length --name  {params.assembly_path} | awk -v  OFS='\t' '{{print $1, $2}}' > {output.length}"


rule make_contig_list_c:
    output:
        contig_txt="output/1_abricate_c/{assembly_name}/{assembly_name}.contig_list"
    params:
        contig_list = lambda wildcards: ASSEMBLY_NAME_CONTIG_LIST_DICT[wildcards.assembly_name]
    run:
        contig_txt=open(output[0],"w")
        for contig in params[0]:
            print(contig, file=contig_txt)
        contig_txt.close()


rule abricate_c:
    input:
        db="output/0_setup_abricate_db/abricate.txt",
        contig_txt="output/1_abricate_c/{assembly_name}/{assembly_name}.contig_list"
    output:
        abricate="output/1_abricate_c/{assembly_name}/{assembly_name}.tsv",
        fasta=temp("output/1_abricate_c/{assembly_name}/{assembly_name}.fasta"),
        length="output/1_abricate_c/{assembly_name}/{assembly_name}.length"
    params:
        assembly_path = lambda wildcards: ASSEMBLY_NAME_PATH_DICT[wildcards.assembly_name]
    conda: "environment.yaml"
    shell:
        "seqkit grep -n -f {input.contig_txt} {params.assembly_path} -o {output.fasta};"
        "./bin/abricate/bin/abricate --db user_db {params.assembly_path} > {output.abricate};"
        "seqkit fx2tab --length --name  {output.fasta} | awk -v  OFS='\t' '{{print $1, $2}}' > {output.length}"


def give_abricate_output(input_format):
    if input_format == "contigs":
        return expand("output/1_abricate_c/{assembly_name}/{assembly_name}.tsv", assembly_name=ASSEMBLY_NAMES)
    elif input_format == "assemblies":
        return expand("output/1_abricate_a/{assembly_name}/{assembly_name}.tsv", assembly_name=ASSEMBLY_NAMES)

def give_length_for_assemblies(input_format):
    if input_format == "contigs":
        return expand("output/1_abricate_c/{assembly_name}/{assembly_name}.length", assembly_name=ASSEMBLY_NAMES)
    elif input_format == "assemblies":
        return expand("output/1_abricate_a/{assembly_name}/{assembly_name}.length", assembly_name=ASSEMBLY_NAMES)


rule abricate_merge:
    input:
        give_abricate_output(config["input_format"])
    output:
         no_head=temp("output/2_filter_gene_observations/1_abricate_all_raw_no_header.tsv"),
         tsv="output/2_filter_gene_observations/1_abricate_all.tsv",
         report="output/2_filter_gene_observations/1_abricate_all.report"
    shell:
        "cat  {input} | grep -v '^#'  > {output.no_head};"
        "cat bin/abricate_header.txt {output.no_head} > {output.tsv};"
        "echo Total number of gene observations in unfiltered abricate results: > {output.report};"
        "cat {output.no_head} | wc -l  >> {output.report}"

rule filter_abricate_quality:
    input:
        tsv="output/2_filter_gene_observations/1_abricate_all.tsv"
    output:
        tsv="output/2_filter_gene_observations/2_abricate_filter_qual.tsv",
        report="output/2_filter_gene_observations/2_abricate_filter_qual.report"
    run:
        # Counter of how many observations that would be discarted at each threshold
        COV_user, COV_100, COV_95, COV_90, COV_85, ID_user, ID_100, ID_95, ID_90, ID_85 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        total_input_observations, total_discarded_observation,  col_index_COV, col_index_ID  = 0, 0, 9, 10 
        user_COV_threshold, user_ID_threshold  = float(config["min_coverage_abricate"]), float(config["min_identity_abricate"])
        
        input_tsv = open("output/2_filter_gene_observations/1_abricate_all.tsv", "r")
        output_tsv = open("output/2_filter_gene_observations/2_abricate_filter_qual.tsv","w")
        for line in input_tsv:
            line = line = line.strip()
            if line.startswith("#"):
                print(line, file = output_tsv)
            else:
                total_input_observations += 1
                line_list = line.split()
                COV_per = float(line_list[col_index_COV])
                ID_per = float(line_list[col_index_ID])
                
                if COV_per < 100:
                    COV_100 +=1
                    if COV_per < 95:
                        COV_95 += 1
                        if COV_per < 90:
                            COV_90 +=1
                            if COV_per < 85: COV_85 += 1
                if COV_per < user_COV_threshold: COV_user += 1
                if ID_per < 100:
                    ID_100 +=1
                    if ID_per < 95:
                        ID_95 += 1
                        if ID_per < 90:
                            ID_90 +=1
                            if ID_per < 85: ID_85 += 1
                if ID_per < user_ID_threshold: ID_user += 1
                if  COV_per >= user_COV_threshold and ID_per >= user_ID_threshold: print(line, file = output_tsv)
                else:
                    total_discarded_observation += 1 
        input_tsv.close
        output_tsv.close
        output_report = open("output/2_filter_gene_observations/2_abricate_filter_qual.report","w")
        a = "Minimum coverage in percentage threshold:" + "\t" + str(user_COV_threshold) + ".\t" + "Sequences excluded due to low coverage:           " + "\t" + str(COV_user)
        b = a +"\nMinimum sequence identity  in percentage threshold::" + "\t" + str(user_ID_threshold )  + ".\t" + "Sequences excluded due to low identity:" + "\t" + str(ID_user)
        c = b + "\n\nNumber of gene observations before quality filtering:" + "\t" + str(total_input_observations)
        d = c + "\nNumber of gene observations after quality filtering:" + "\t" + str(total_input_observations - total_discarded_observation)
        e = d + "\nNumber of gene observations discarded by quality filtering:" + "\t" + str(total_discarded_observation)
        f = e +"\n\n\nHypothetical number of sequences discarded due to low coverage if the following thresholds were used"
        g = f + "\nCOV_100" + "\t" + "COV_95" + "\t" + "COV_90" + "\t"  + "COV_85\n" +  str(COV_100) + "\t" + str(COV_95) + "\t" + str(COV_90) + "\t"  + str(COV_85)
        h = g + "\n\nHypothetical number of sequences discarded due to low sequence identity if the following thresholds were used"
        i = h + "\nID_100" + "\t" + "ID_95" + "\t" + "ID_90" + "\t"  + "ID_85\n" +  str(ID_100) + "\t" + str(ID_95) + "\t" + str(ID_90) + "\t"  + str(ID_85)
        print(i, file = output_report)
        output_report.close

rule merge_length_files:
    input:
        give_length_for_assemblies(config["input_format"])
    output:
        temp("output/2_filter_gene_observations/3_all_contig.lengths")
    shell:
        "cat {input} > {output}"


rule filter_abricate_space_for_flanks:
    input:
        length="output/2_filter_gene_observations/3_all_contig.lengths",
        tsv="output/2_filter_gene_observations/2_abricate_filter_qual.tsv"
    output:
        tsv="output/2_filter_gene_observations/3_final_gene_results.tsv",
        report="output/2_filter_gene_observations/3_abricate_filter_length.report"
    run:
        LENGTH_DICT = {}
        counter_accepted = 0
        input_length = open("output/2_filter_gene_observations/3_all_contig.lengths", "r")
        for line in input_length:
            line = line.strip()
            line_list = line.split()
            LENGTH_DICT[line_list[0]] = int(line_list[1])
        input_length.close
        
        # Counter of how many observations that would be discarted at each threshold
        UP_user, UP_half, UP_double, DOWN_user, DOWN_half, DOWN_double = 0, 0, 0, 0, 0, 0 
        total_input_observations, total_discarded_observation = 0, 0  
        col_index_ASSEMBLY_PATH, col_index_SEQ, col_index_START, col_index_END, col_index_STRAND  = 0, 1, 2, 3, 4
        user_upstreams=float(config["flank_length_upstreams"])
        user_downstreams=float(config["flank_length_downstreams"])         

        input_tsv = open("output/2_filter_gene_observations/2_abricate_filter_qual.tsv", "r")
        output_tsv = open("output/2_filter_gene_observations/3_final_gene_results.tsv", "w")
        for line in input_tsv:
            line = line.strip()
            if line.startswith("#"): print(line, file = output_tsv)
            else:
                space_for_flank_up, space_for_flank_down = 0, 0
                total_input_observations += 1
                line_list = line.split()
                path = line_list[col_index_ASSEMBLY_PATH]
                seq_name = line_list[col_index_SEQ]
                gene_start = int(line_list[col_index_START])
                gene_end = int(line_list[col_index_END])
                gene_strand = line_list[col_index_STRAND]
                if gene_strand == "+":
                    space_for_flank_up = gene_start - 1
                    space_for_flank_down = int(LENGTH_DICT[seq_name]) - gene_end
                if gene_strand == "-":
                    space_for_flank_down = gene_start - 1
                    space_for_flank_up = int(LENGTH_DICT[seq_name]) - gene_end
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
                    assembly_name = PATH_ASSEMBLY_NAME_DICT[path]
                    new_line = line + "\t" + assembly_name + "\t" + "i" + str(counter_accepted)
                    print(new_line, file = output_tsv)
                else:
                    total_discarded_observation += 1
        input_tsv.close
        output_tsv.close
        output_report = open("output/2_filter_gene_observations/3_abricate_filter_length.report", "w")
        a = "Flank length upstreams:" + "\t" + str(int(user_upstreams))
        b = a +"\nFlank length downstreams:" + "\t" + str(int(user_downstreams))
        c = b + "\n\nNumber of gene observations before flank space filtering:" + "\t" + str(total_input_observations)
        d = c + "\nNumber of gene observations after flank space filtering:" + "\t" + str(total_input_observations - total_discarded_observation)
        e = d + "\nNumber of gene observations discarded by flank space filtering:" + "\t" + str(total_discarded_observation)
        f = e +"\n\n\nHypothetical number of sequences discarded due to insufficient upstreams flank space if the following thresholds were used"
        g = f + "\n" + str(round(user_upstreams / 2)) + "\t" + str(round(user_upstreams)) + "\t"  + str(round(user_upstreams * 2)) 
        h = g + "\n" +  str(UP_half) + "\t" + str(UP_user) + "\t"  + str(UP_double)
        i = h + "\n\nHypothetical number of sequences discarded due to insufficient downstreams flank space if the following thresholds were used"
        j = i + "\n" + str(round(user_downstreams / 2)) + "\t" + str(round(user_downstreams)) + "\t"  + str(round(user_downstreams* 2))
        k = j + "\n" +  str(DOWN_half) + "\t" + str(DOWN_user) + "\t"  + str(DOWN_double)
        print(k, file = output_report)
        output_report.close()



rule add_file_name_to_tsv:
    input:
        "output/2_filter_gene_observations/3_final_gene_results.tsv"
    output:
        temp("output/3_extract_sequences/final_gene_results_with_future_filename.tsv")
    run:
        input=open(input[0], "r")
        output=open(output[0], "w")
        for line in input:
            line = line.strip()
            line_list = line.split()
            ASSEMBLY_NAME = line_list[13]
            future_file_name = ASSEMBLY_NAME + ".tsv"
            new_line = line + "\t" + future_file_name
            print(new_line, file=output) 
        input.close()
        output.close()
        
checkpoint split_abricate_results:
    input:
        "output/3_extract_sequences/final_gene_results_with_future_filename.tsv"
    output:
        clusters=directory("output/3_extract_sequences/1_abricate_results_per_assembly")
    shell:
        "mkdir output/3_extract_sequences/1_abricate_results_per_assembly;"
        "cd output/3_extract_sequences/1_abricate_results_per_assembly;"
        "cat ../final_gene_results_with_future_filename.tsv | grep -v '^#' | awk  '{{print>$16}}'"
# split into files named as assembly_name



rule make_bedfiles:
    input:
        "output/3_extract_sequences/1_abricate_results_per_assembly/{a}.tsv"
    output:
        "output/3_extract_sequences/2_bedfiles/flanks_with_gene/{a}/{a}.bed",
        "output/3_extract_sequences/2_bedfiles/just_gene/{a}/{a}.bed",
        "output/3_extract_sequences/2_bedfiles/masked_gene/{a}/{a}.bed"
    run:
        input=open(input[0], "r")
        output_flanks_with_gene=open(output[0], "w")
        output_just_gene=open(output[1], "w")
        output_masked_gene=open(output[2], "w")
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

            # for flanks_with_gene - whole region - with bedtools getfasta
            if gene_strand == "+":
                START_flanks_with_gene = START_old - user_upstreams -1
                END_flanks_with_gene = END_old + user_downstreams
            if gene_strand == "-":
                START_flanks_with_gene = START_old - user_downstreams -1
                END_flanks_with_gene = END_old + user_upstreams
            out_line = contig_name +  "\t" + str(START_flanks_with_gene) +  "\t" + str(END_flanks_with_gene)
            out_line = out_line + "\t" + ID + "\t" + "1" + "\t" + gene_strand
            print(out_line, file=output_flanks_with_gene)
 
            #just_gene - just cut out the gene  with bedtools getfasta
            out_line = contig_name +  "\t" + str(START_old - 1) +  "\t" + str(END_old)
            out_line = out_line + "\t" + ID + "\t" + "1" + "\t" + gene_strand
            print(out_line, file=output_just_gene)

            #masked_gene - made with maskfasta from flanks_with_gene_fasta. All strands will be positive
            START_masked_gene = user_upstreams + 1 - 1
            gene_length = END_old - START_old + 1
            END_masked_gene = user_upstreams + gene_length
            out_line = ID +  "\t" + str(START_masked_gene) + "\t" + str(END_masked_gene)
            print(out_line, file =output_masked_gene )
        input.close()
        output_flanks_with_gene.close
        output_just_gene.close
        output_masked_gene.close
        


rule bedtools_flanks_with_gene:
    input:
        bed="output/3_extract_sequences/2_bedfiles/flanks_with_gene/{a}/{a}.bed"
    output:
        fasta="output/3_extract_sequences/3_gene_fasta_files_per_assembly/flanks_with_gene/{a}/{a}.fa"
    conda: "environment.yaml"
    params:
        assembly_path=lambda wildcards: ASSEMBLY_NAME_PATH_DICT[wildcards.a]
    shell:
        '''
        bedtools getfasta -fi {params.assembly_path} -bed {input.bed} -s -nameOnly > {output.fasta};
        sed -i -e "s/[(|)|+]//g" -e "s/-//g" {output.fasta}
        '''

rule bedtools_just_gene:
    input:
        whole_flank="output/3_extract_sequences/4_gene_fasta_files_pooled/flanks_with_gene.fa",
        bed="output/3_extract_sequences/2_bedfiles/just_gene/{a}/{a}.bed"
    output:
        fasta="output/3_extract_sequences/3_gene_fasta_files_per_assembly/just_gene/{a}/{a}.fa"
    conda: "environment.yaml"
    params:
        assembly_path=lambda wildcards: ASSEMBLY_NAME_PATH_DICT[wildcards.a]
    shell:
        '''
        bedtools getfasta -fi {params.assembly_path} -bed {input.bed} -s -nameOnly > {output.fasta};
        sed -i -e "s/[(|)|+]//g" -e "s/-//g" {output.fasta}
        ''' 


rule bedtools_masked_gene:
    input:
        just_gene="output/3_extract_sequences/4_gene_fasta_files_pooled/just_gene.fa",
        bed="output/3_extract_sequences/2_bedfiles/masked_gene/{a}/{a}.bed",
        fasta="output/3_extract_sequences/3_gene_fasta_files_per_assembly/flanks_with_gene/{a}/{a}.fa"
    output:
        fasta="output/3_extract_sequences/3_gene_fasta_files_per_assembly/masked_gene/{a}/{a}.fa"
    conda: "environment.yaml"
    shell:
        "bedtools maskfasta -fi {input.fasta} -bed {input.bed} -fo {output.fasta}"



def flanks_with_gene_input(wildcards):
    checkpoint_output = checkpoints.split_abricate_results.get(**wildcards).output[0]
    return expand("output/3_extract_sequences/3_gene_fasta_files_per_assembly/flanks_with_gene/{a}/{a}.fa",
           a=glob_wildcards(os.path.join(checkpoint_output, "{a}.tsv")).a)



def just_gene_input(wildcards):
    checkpoint_output = checkpoints.split_abricate_results.get(**wildcards).output[0]
    return expand("output/3_extract_sequences/3_gene_fasta_files_per_assembly/just_gene/{a}/{a}.fa",
           a=glob_wildcards(os.path.join(checkpoint_output, "{a}.tsv")).a)


def masked_gene_input(wildcards):
    checkpoint_output = checkpoints.split_abricate_results.get(**wildcards).output[0]
    return expand("output/3_extract_sequences/3_gene_fasta_files_per_assembly/masked_gene/{a}/{a}.fa",
           a=glob_wildcards(os.path.join(checkpoint_output, "{a}.tsv")).a)


rule pool_flanks_with_gene:
    input:
        flanks_with_gene_input
    output:
        temp=temp("output/3_extract_sequences/4_gene_fasta_files_pooled/temp_flanks_with_gene.fa"),
        fasta="output/3_extract_sequences/4_gene_fasta_files_pooled/flanks_with_gene.fa"
    conda: "environment.yaml"
    shell:
        "cat {input}  > {output.temp};"
        "cat {output.temp} | seqkit sort -N -i -2  > {output.fasta}"


rule pool_just_gene:
    input:
        just_gene_input
    output:
        temp=temp("output/3_extract_sequences/4_gene_fasta_files_pooled/temp_just_gene.fa"),
        fasta="output/3_extract_sequences/4_gene_fasta_files_pooled/just_gene.fa"
    conda: "environment.yaml"
    shell:
        "cat {input}  > {output.temp};"
        "cat {output.temp} | seqkit sort -N -i -2  > {output.fasta}"


rule pool_masked_gene:
    input:
        masked_gene_input
    output:
        temp=temp("output/3_extract_sequences/4_gene_fasta_files_pooled/temp_masked_gene.fa"),
        fasta="output/3_extract_sequences/4_gene_fasta_files_pooled/masked_gene.fa"
    conda: "environment.yaml"
    shell:
        "cat {input}  > {output.temp};"
        "cat {output.temp} | seqkit sort -N -i -2  > {output.fasta}"


rule extract_relavant_reference_sequences:
    input:
        flanks="output/3_extract_sequences/4_gene_fasta_files_pooled/flanks_with_gene.fa",
        just_gene="output/3_extract_sequences/4_gene_fasta_files_pooled/just_gene.fa",
        masked="output/3_extract_sequences/4_gene_fasta_files_pooled/masked_gene.fa",
        ref_seq="bin/abricate/db/user_db/sequences",
        final_results="output/2_filter_gene_observations/3_final_gene_results.tsv"
    output:
        list=temp("output/4_cluster_by_gene_family/final_output_list.txt"),
        fasta="output/4_cluster_by_gene_family/relevant_ref_seqs.fasta"
    conda: "environment.yaml"
    shell:
        "awk '{{print $13}}' {input.final_results} |awk 'NR>1' | sort -u > {output.list};"
        "seqkit grep -n -f {output.list} {input.ref_seq} > {output.fasta}"


rule cd_hit:
    input:
        "output/4_cluster_by_gene_family/relevant_ref_seqs.fasta"
    output:
        clus="output/4_cluster_by_gene_family/cd_hit.clstr",
        genes=temp("output/4_cluster_by_gene_family/cd_hit"),
        new_head=temp("output/4_cluster_by_gene_family/cd_hit.txt"),
        temp=temp("output/4_cluster_by_gene_family/cd_hit.temp")
    conda: "environment.yaml"
    params:
        c= config["cluster_identity_cd_hit"],
        n= config["cluster_wordsize_cd_hit"],
        s= config["cluster_length_dif_cd_hit"]
    shell:
        '''
        cd-hit-est -i {input} -o {output.genes} -M 0 -d 0 -c {params.c} -n {params.n} -s {params.s}  -sc 1 -g 1;
        perl bin/clstr2txt.pl {output.clus} | tail -n+2 > {output.temp};
        cat bin/cd-hit-txt_header.txt {output.temp} > {output.new_head}
        '''


rule process_cd_hit_results:
    input: 
        "output/4_cluster_by_gene_family/cd_hit.txt"
    output:
        "output/4_cluster_by_gene_family/cd_hit.processed"
    run:
        dict = {}
        input=open(input[0], "r")
        output=open(output[0], "w")

        for line in input:
            if not line.startswith("#"):
                line = line.strip()
                line_list = line.split()
                gene_name = line_list[0]
                clus_number = line_list[1]
                if clus_number not in dict:
                    # Split name so only first part of name is used
                    first_part_gene_name = re.split(', |\'|\/|\(|\)|_| ', gene_name)[0]
                    dict[clus_number] = clus_number + "_" + first_part_gene_name + ".tsv"
                new_line = line + "\t" + dict[clus_number]
                print(new_line, file=output)        
        input.close()
        output.close()


checkpoint split_cd_hit_results:
    input:
        "output/4_cluster_by_gene_family/cd_hit.processed"
    output:
        clusters=directory("output/4_cluster_by_gene_family/cd_hit_per_cluster")
    shell:
        "mkdir output/4_cluster_by_gene_family/cd_hit_per_cluster;"
        "cd output/4_cluster_by_gene_family/cd_hit_per_cluster;"
        "cat ../cd_hit.processed | awk  '{{print>$8}}'"

## 5 ########################################################


rule extract_cluster_rows_from_results_file:
    input:
        results="output/2_filter_gene_observations/3_final_gene_results.tsv",
        list="output/4_cluster_by_gene_family/cd_hit_per_cluster/{c}.tsv"
    output:
        results="output/5_gene_clusters/{c}/{c}.tsv",
        temp=temp("output/5_gene_clusters/{c}/{c}.temp")
    params:
        awk="'NR==FNR{a[$1];next}($13 in a){{print $0}}'"
    shell:
        '''
        awk {params.awk} {input.list} {input.results} > {output.temp};
        cat bin/abricate_header.txt {output.temp} > {output.results}
        '''


rule extract_cluster_fastas:
    input:
        tsv="output/5_gene_clusters/{c}/{c}.tsv",
        flanks_with_gene_fasta="output/3_extract_sequences/4_gene_fasta_files_pooled/flanks_with_gene.fa",
        just_gene_fasta="output/3_extract_sequences/4_gene_fasta_files_pooled/just_gene.fa",
        masked_gene_fasta="output/3_extract_sequences/4_gene_fasta_files_pooled/masked_gene.fa"
    output:
        ID_list=temp("output/5_gene_clusters/{c}/{c}.ID_list"),
        flanks_with_gene_fasta="output/5_gene_clusters/{c}/{c}.flanks_with_gene_fasta",
        just_gene_fasta="output/5_gene_clusters/{c}/{c}.just_gene_fasta",
        masked_gene_fasta="output/5_gene_clusters/{c}/{c}.masked_gene_fasta"
    conda: "environment.yaml"
    shell:
        "awk '{{print $15}}' {input.tsv} > {output.ID_list};"
        "seqkit grep -n -f {output.ID_list} {input.flanks_with_gene_fasta} > {output.flanks_with_gene_fasta};"
        "seqkit grep -n -f {output.ID_list} {input.just_gene_fasta} > {output.just_gene_fasta};"
        "seqkit grep -n -f {output.ID_list} {input.masked_gene_fasta} > {output.masked_gene_fasta}"


  

def split_just_gene_fasta_file_input(wildcards):
    checkpoint_output = checkpoints.split_cd_hit_results.get(**wildcards).output[0]
    return expand("output/5_gene_clusters/{c}/{c}.just_gene_fasta",
           c=glob_wildcards(os.path.join(checkpoint_output, "{c}.tsv")).c)


rule translate_reference_database:
    input:
        split_just_gene_fasta_file_input
    output:
        temp("output/user_db")
    conda: "environment.yaml"
    params:
        fasta=config["database"]
    shell:
        "seqkit translate {params.fasta} > {output}"


rule prokka:
    input:
        AA_db="output/user_db",
        fasta="output/5_gene_clusters/{c}/{c}.flanks_with_gene_fasta"
    output:
        gff="output/5_gene_clusters/{c}/prokka/{c}.gff"
    params:
        outdir=lambda wildcards: "output/5_gene_clusters/" +  wildcards.c + "/prokka",
        c=lambda wildcards: wildcards.c
    conda: "environment.yaml" 
    shell:
        "prokka --proteins {input.AA_db} --outdir {params.outdir} --force --prefix {params.c} {input.fasta}"



rule make_gggenes_input:
    input:
        "output/5_gene_clusters/{c}/prokka/{c}.gff"
    output:
        "output/5_gene_clusters/{c}/{c}.gggenes"
    shell:
        "python3 ./bin/gff_to_gggenes.py {input} {output}"



rule kma_index:
    input:
        gggenes="output/5_gene_clusters/{c}/{c}.gggenes",
        masked="output/5_gene_clusters/{c}/{c}.masked_gene_fasta",
        just_gene="output/5_gene_clusters/{c}/{c}.just_gene_fasta"
    output:
        "output/5_gene_clusters/{c}/kma_index/masked_gene.comp.b",
        "output/5_gene_clusters/{c}/kma_index/masked_gene.length.b",
        "output/5_gene_clusters/{c}/kma_index/masked_gene.name",
        "output/5_gene_clusters/{c}/kma_index/masked_gene.seq.b",
        "output/5_gene_clusters/{c}/kma_index/just_gene.comp.b",
        "output/5_gene_clusters/{c}/kma_index/just_gene.length.b",
        "output/5_gene_clusters/{c}/kma_index/just_gene.name",
        "output/5_gene_clusters/{c}/kma_index/just_gene.seq.b"
    conda: "environment.yaml"
    params:
        outfolder_masked=lambda wildcards: "output/5_gene_clusters/" +  wildcards.c + "/kma_index/masked_gene",
        outfolder_just_gene=lambda wildcards: "output/5_gene_clusters/" +  wildcards.c + "/kma_index/just_gene"
    shell:
        "kma index -i {input.masked} -o {params.outfolder_masked};"
        "kma index -i {input.just_gene} -o {params.outfolder_just_gene}"




rule kma_dist:
    input:
        m="output/5_gene_clusters/{c}/kma_index/masked_gene.seq.b",
        j="output/5_gene_clusters/{c}/kma_index/just_gene.seq.b"
    output:
        m="output/5_gene_clusters/{c}/{c}.masked_gene_dist",
        j="output/5_gene_clusters/{c}/{c}.just_gene_dist"
    params:
        m_db=lambda wildcards: "output/5_gene_clusters/" +  wildcards.c + "/kma_index/masked_gene",
        j_db=lambda wildcards: "output/5_gene_clusters/" +  wildcards.c + "/kma_index/just_gene",
        m_dist_measure=config["distance_measure_flanks_masked_gene"],
        j_dist_measure=config["distance_measure_just_gene"]
    conda: "environment.yaml"
    shell:
        "kma dist -t_db {params.m_db} -o {output.m} -d {params.m_dist_measure};"
        "kma dist -t_db {params.j_db} -o {output.j} -d {params.j_dist_measure}"



rule clearcut:
    input:
        m="output/5_gene_clusters/{c}/{c}.masked_gene_dist",
        j="output/5_gene_clusters/{c}/{c}.just_gene_dist"
    output:
        m="output/5_gene_clusters/{c}/{c}.masked_gene_tree",
        j="output/5_gene_clusters/{c}/{c}.just_gene_tree"

    conda: "environment.yaml"
    shell:
        '''
        echo To few leaves to make a tree > {output.m};
        echo To few leaves to make a tree > {output.j};
        clearcut --in={input.m} --out={output.m} --neighbor || true ;
        clearcut --in={input.j} --out={output.j} --neighbor || true 
        '''
        


def tree_input(wildcards):
    checkpoint_output = checkpoints.split_cd_hit_results.get(**wildcards).output[0]
    return expand("output/5_gene_clusters/{c}/{c}.masked_gene_tree",
           c=glob_wildcards(os.path.join(checkpoint_output, "{c}.tsv")).c)


## 99 #############################################


rule get_trees:
    input:
        tree_input
    output:
        "output/99_trees_and_distance_matrixes_done"
    shell:
        "echo done > {output}"




rule finito:
    output:
        "output/99_config.yaml"
    shell:
        "cp config.yaml {output}"

