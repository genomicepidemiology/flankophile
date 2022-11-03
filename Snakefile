# FLANKOPHILE version 0.1.6
# Alex Vincent Thorn

configfile: "config.yaml"

import os
import re

ASSEMBLY_NAMES = []
ASSEMBLY_NAME_PATH_DICT = {}
PATH_ASSEMBLY_NAME_DICT = {}

with open(config["input_list"], 'r') as file:
    for line in file:
        line = line.strip()
        if not line.startswith("#"):
            line_list = line.split("\t")
            assembly_name = line_list[0]
            assembly_path = line_list[1]
            if assembly_name not in ASSEMBLY_NAMES:
                ASSEMBLY_NAMES.append(assembly_name)
                ASSEMBLY_NAME_PATH_DICT[assembly_name] = assembly_path
                PATH_ASSEMBLY_NAME_DICT[assembly_path] = assembly_name

            
rule all:
    input: 
       "output/2_filter/3_final_gene_results.tsv",
       "output/3_define_clusters/cd_hit.processed",
       "output/99_config.yaml"


####################################################################################

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
         tsv="output/1_search/search_raw.tsv"
    shell:
        "cat  {input} | grep -v '^#'| sed 's/\t\t/\t/g'   > {output.no_head};"
        "cat bin/abricate_header.txt {output.no_head} > {output.tsv};"


         

rule filter_abricate_quality:
    input:
        tsv="output/1_search/search_raw.tsv"
    output:
        tsv="output/2_filter/2_abricate_filter_qual.tsv"
    run:
        input_tsv = open("output/1_search/search_raw.tsv", "r")
        output_tsv = open("output/2_filter/2_abricate_filter_qual.tsv","w")
        row_count = 0
        for line in input_tsv:
            line = line = line.strip()
            if line.startswith("#"):
                print(line, file = output_tsv)
            else:
                row_count += 1
                line_list = line.split("\t")
                assembly_name = PATH_ASSEMBLY_NAME_DICT[line_list[0]]


                new_line = line + "\t" + assembly_name + "\t" + "i" + str(row_count)
                print(new_line, file = output_tsv)

        input_tsv.close
        output_tsv.close
        




rule filter_abricate_space_for_flanks:
    input:
        tsv="output/2_filter/2_abricate_filter_qual.tsv"
    output:
        tsv="output/2_filter/3_final_gene_results.tsv",
        report="output/2_filter/3_abricate_filter_length.report"
    run:
        counter_accepted = 0

        # Counter of how many observations that would be discarted at each threshold
        UP_user, UP_half, UP_double, DOWN_user, DOWN_half, DOWN_double = 0, 0, 0, 0, 0, 0 
        
        total_input_observations, total_discarded_observation = 0, 0  
        col_index_ASSEMBLY_PATH, col_index_SEQ, col_index_START, col_index_END, col_index_STRAND  = 0, 1, 2, 3, 4
        user_upstreams=float(config["flank_length_upstreams"])
        user_downstreams=float(config["flank_length_downstreams"])
                 

        input_tsv = open("output/2_filter/2_abricate_filter_qual.tsv", "r")
        output_tsv = open("output/2_filter/3_final_gene_results.tsv", "w")
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
                
                assembly_name = PATH_ASSEMBLY_NAME_DICT[path]
                
                path_to_length_file = "output/1_search/contig_length/" + assembly_name + "/" + assembly_name + ".length"
                
                this_contig_len = "unknown"
                        
                input_length = open(path_to_length_file, "r")
                for line_l in input_length:
                    if this_contig_len == "unknown":
                        line_l = line_l.strip()
                        line_l_list = line_l.split("\t")
                        if line_l_list[0] == seq_name:
                            this_contig_len = int(line_l_list[1])
                input_length.close  
                    
                
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
                   # new_line = line + "\t" + assembly_name + "\t" + "i" + str(counter_accepted)
                    #new_line = line + "\t" + "i" + str(counter_accepted)
                    #print(new_line, file = output_tsv)
                    print(line, file = output_tsv)
                else:
                    total_discarded_observation += 1
        input_tsv.close
        output_tsv.close
        output_report = open("output/2_filter/3_abricate_filter_length.report", "w")
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
        "output/2_filter/3_final_gene_results.tsv"
    output:
        temp("output/temp_1/final_gene_results_with_future_filename.tsv")
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
        "output/temp_1/final_gene_results_with_future_filename.tsv"
    output:
        clusters=directory("output/temp_1/1_abricate_results_per_assembly")
    shell:
        "mkdir output/temp_1/1_abricate_results_per_assembly;"
        "cd output/temp_1/1_abricate_results_per_assembly;"
        "cat ../final_gene_results_with_future_filename.tsv | grep -v '^#' | awk  '{{print>$16}}'"


rule make_bedfiles:
    input:
        "output/temp_1/1_abricate_results_per_assembly/{a}.tsv"
    output:
        temp("output/temp_2/bedfiles/flanks_with_gene/{a}/{a}.bed"),
        temp("output/temp_2/bedfiles/just_gene/{a}/{a}.bed"),
        temp("output/temp_2/bedfiles/masked_gene/{a}/{a}.bed")
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
        bed="output/temp_2/bedfiles/flanks_with_gene/{a}/{a}.bed"
    output:
        fasta=temp("output/temp_2/gene_fasta_files_per_assembly/flanks_with_gene/{a}/{a}.fa")
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
        whole_flank="output/2_filter/4_gene_fasta_files_pooled/flanks_with_gene.fa",
        bed="output/temp_2/bedfiles/just_gene/{a}/{a}.bed"
    output:
        fasta=temp("output/temp_2/gene_fasta_files_per_assembly/just_gene/{a}/{a}.fa")
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
        just_gene="output/2_filter/4_gene_fasta_files_pooled/just_gene.fa",
        bed="output/temp_2/bedfiles/masked_gene/{a}/{a}.bed",
        fasta="output/temp_2/gene_fasta_files_per_assembly/flanks_with_gene/{a}/{a}.fa"
    output:
        fasta=temp("output/temp_2/gene_fasta_files_per_assembly/masked_gene/{a}/{a}.fa")
    conda: "environment.yaml"
    shell:
        "bedtools maskfasta -fi {input.fasta} -bed {input.bed} -fo {output.fasta}"



def flanks_with_gene_input(wildcards):
    checkpoint_output = checkpoints.split_abricate_results.get(**wildcards).output[0]
    return expand("output/temp_2/gene_fasta_files_per_assembly/flanks_with_gene/{a}/{a}.fa",
           a=glob_wildcards(os.path.join(checkpoint_output, "{a}.tsv")).a)



def just_gene_input(wildcards):
    checkpoint_output = checkpoints.split_abricate_results.get(**wildcards).output[0]
    return expand("output/temp_2/gene_fasta_files_per_assembly/just_gene/{a}/{a}.fa",
           a=glob_wildcards(os.path.join(checkpoint_output, "{a}.tsv")).a)


def masked_gene_input(wildcards):
    checkpoint_output = checkpoints.split_abricate_results.get(**wildcards).output[0]
    return expand("output/temp_2/gene_fasta_files_per_assembly/masked_gene/{a}/{a}.fa",
           a=glob_wildcards(os.path.join(checkpoint_output, "{a}.tsv")).a)


rule pool_flanks_with_gene:
    input:
        flanks_with_gene_input
    output:
        temp=temp("output/temp_2/gene_fasta_files_pooled/temp_flanks_with_gene.fa"),
        fasta="output/2_filter/4_gene_fasta_files_pooled/flanks_with_gene.fa"
    conda: "environment.yaml"
    shell:
        "cat {input}  > {output.temp};"
        "cat {output.temp} | seqkit sort -N -i -2  > {output.fasta};"
        "rm -r output/temp_1"


rule pool_just_gene:
    input:
        just_gene_input
    output:
        temp=temp("output/temp_2/gene_fasta_files_pooled/temp_just_gene.fa"),
        fasta="output/2_filter/4_gene_fasta_files_pooled/just_gene.fa"
    conda: "environment.yaml"
    shell:
        "cat {input}  > {output.temp};"
        "cat {output.temp} | seqkit sort -N -i -2  > {output.fasta}"


rule pool_masked_gene:
    input:
        masked_gene_input
    output:
        temp=temp("output/temp_2/gene_fasta_files_pooled/temp_masked_gene.fa"),
        fasta="output/2_filter/4_gene_fasta_files_pooled/masked_gene.fa"
    conda: "environment.yaml"
    shell:
        "cat {input}  > {output.temp};"
        "cat {output.temp} | seqkit sort -N -i -2  > {output.fasta}"


rule extract_relavant_reference_sequences:
    input:
        flanks="output/2_filter/4_gene_fasta_files_pooled/flanks_with_gene.fa",
        just_gene="output/2_filter/4_gene_fasta_files_pooled/just_gene.fa",
        masked="output/2_filter/4_gene_fasta_files_pooled/masked_gene.fa",
        final_results="output/2_filter/3_final_gene_results.tsv"
    output:
        list=temp("output/3_define_clusters/final_output_list.txt"),
        fasta=temp("output/3_define_clusters/relevant_ref_seqs.fasta")
    conda: "environment.yaml"
    shell:
        "awk '{{print $13}}' {input.final_results} |awk 'NR>1' | sort -u > {output.list};"
        "seqkit grep -n -f {output.list} bin/abricate/db/user_db/sequences > {output.fasta}"


rule cd_hit:
    input:
        "output/3_define_clusters/relevant_ref_seqs.fasta"
    output:
        clus="output/3_define_clusters/cd_hit.clstr",
        genes=temp("output/3_define_clusters/cd_hit"),
        new_head=temp("output/3_define_clusters/cd_hit.txt"),
        temp=temp("output/3_define_clusters/cd_hit.temp")
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
        "output/3_define_clusters/cd_hit.txt"
    output:
        "output/3_define_clusters/cd_hit.processed"
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
                    # first_part_gene_name = re.split(', |\'|\/|\(|\)|_| ', gene_name)[0]
                    first_part_gene_name = re.sub(r'\W', '_', gene_name)
                    dict[clus_number] = clus_number + "_" + first_part_gene_name + ".tsv"
                new_line = line + "\t" + dict[clus_number]
                print(new_line, file=output)        
        input.close()
        output.close()


checkpoint split_cd_hit_results:
    input:
        "output/3_define_clusters/cd_hit.processed"
    output:
        clusters=directory("output/3_define_clusters/cd_hit_per_cluster")
    shell:
        "mkdir output/3_define_clusters/cd_hit_per_cluster;"
        "cd output/3_define_clusters/cd_hit_per_cluster;"
        "cat ../cd_hit.processed | awk  '{{print>$8}}'"




rule extract_cluster_rows_from_results_file:
    input:
        results="output/2_filter/3_final_gene_results.tsv",
        list="output/3_define_clusters/cd_hit_per_cluster/{c}.tsv"
    output:
        results="output/4_cluster_results/{c}/{c}.tsv",
        temp=temp("output/4_cluster_results/{c}/{c}.temp")
    params:
        awk="'NR==FNR{a[$1];next}($13 in a){{print $0}}'"
    shell:
        '''
        awk {params.awk} {input.list} {input.results} > {output.temp};
        cat bin/abricate_header.txt {output.temp} > {output.results}
        '''


rule extract_cluster_fastas:
    input:
        tsv="output/4_cluster_results/{c}/{c}.tsv",
        flanks_with_gene_fasta="output/2_filter/4_gene_fasta_files_pooled/flanks_with_gene.fa",
        just_gene_fasta="output/2_filter/4_gene_fasta_files_pooled/just_gene.fa",
        masked_gene_fasta="output/2_filter/4_gene_fasta_files_pooled/masked_gene.fa"
    output:
        ID_list=temp("output/4_cluster_results/{c}/{c}.ID_list"),
        flanks_with_gene_fasta="output/4_cluster_results/{c}/{c}.flanks_with_gene_fasta",
        just_gene_fasta="output/4_cluster_results/{c}/{c}.just_gene_fasta",
        masked_gene_fasta="output/4_cluster_results/{c}/{c}.masked_gene_fasta"
    conda: "environment.yaml"
    shell:
        "awk '{{print $15}}' {input.tsv} > {output.ID_list};"
        "seqkit grep -n -f {output.ID_list} {input.flanks_with_gene_fasta} > {output.flanks_with_gene_fasta};"
        "seqkit grep -n -f {output.ID_list} {input.just_gene_fasta} > {output.just_gene_fasta};"
        "seqkit grep -n -f {output.ID_list} {input.masked_gene_fasta} > {output.masked_gene_fasta}"


  

def split_just_gene_fasta_file_input(wildcards):
    checkpoint_output = checkpoints.split_cd_hit_results.get(**wildcards).output[0]
    return expand("output/4_cluster_results/{c}/{c}.just_gene_fasta",
           c=glob_wildcards(os.path.join(checkpoint_output, "{c}.tsv")).c)


rule translate_reference_database:
    input:
        split_just_gene_fasta_file_input
    output:
        temp("output/user_db")
    conda: "environment.yaml"
    shell:
        "seqkit translate bin/abricate/db/user_db/sequences > {output}"


rule prokka:
    input:
        AA_db="output/user_db",
        fasta="output/4_cluster_results/{c}/{c}.flanks_with_gene_fasta"
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
        "output/4_cluster_results/{c}/prokka/{c}.gff"
    output:
        "output/4_cluster_results/{c}/{c}.gggenes"
    shell:
        "python3 ./bin/gff_to_gggenes.py {input} {output}"



rule kma_index:
    input:
        gggenes="output/4_cluster_results/{c}/{c}.gggenes",
        flanks_with_gene="output/4_cluster_results/{c}/{c}.flanks_with_gene_fasta",
        masked="output/4_cluster_results/{c}/{c}.masked_gene_fasta",
        just_gene="output/4_cluster_results/{c}/{c}.just_gene_fasta"
    output:
        temp("output/4_cluster_results/{c}/kma_index/flanks_with_gene.comp.b"),
        temp("output/4_cluster_results/{c}/kma_index/flanks_with_gene.length.b"),
        temp("output/4_cluster_results/{c}/kma_index/flanks_with_gene.name"),
        temp("output/4_cluster_results/{c}/kma_index/flanks_with_gene.seq.b"),
        temp("output/4_cluster_results/{c}/kma_index/masked_gene.comp.b"),
        temp("output/4_cluster_results/{c}/kma_index/masked_gene.length.b"),
        temp("output/4_cluster_results/{c}/kma_index/masked_gene.name"),
        temp("output/4_cluster_results/{c}/kma_index/masked_gene.seq.b"),
        temp("output/4_cluster_results/{c}/kma_index/just_gene.comp.b"),
        temp("output/4_cluster_results/{c}/kma_index/just_gene.length.b"),
        temp("output/4_cluster_results/{c}/kma_index/just_gene.name"),
        temp("output/4_cluster_results/{c}/kma_index/just_gene.seq.b")
    conda: "environment.yaml"
    params:
        outfolder_flanks_with_gene=lambda wildcards: "output/4_cluster_results/" +  wildcards.c + "/kma_index/flanks_with_gene",
        outfolder_masked=lambda wildcards: "output/4_cluster_results/" +  wildcards.c + "/kma_index/masked_gene",
        outfolder_just_gene=lambda wildcards: "output/4_cluster_results/" +  wildcards.c + "/kma_index/just_gene",
        k=config["Kmersize_kma"]
    shell:
        "kma index -k {params.k} -i {input.flanks_with_gene} -o {params.outfolder_flanks_with_gene};"
        "kma index -k {params.k} -i {input.masked} -o {params.outfolder_masked};"
        "kma index -k {params.k} -i {input.just_gene} -o {params.outfolder_just_gene}"




rule kma_dist:
    input:
        ja="output/4_cluster_results/{c}/kma_index/just_gene.comp.b",
        jb="output/4_cluster_results/{c}/kma_index/just_gene.length.b",
        jc="output/4_cluster_results/{c}/kma_index/just_gene.name",
        fa="output/4_cluster_results/{c}/kma_index/flanks_with_gene.comp.b",
        fb="output/4_cluster_results/{c}/kma_index/flanks_with_gene.length.b",
        fc="output/4_cluster_results/{c}/kma_index/flanks_with_gene.name",
        ma="output/4_cluster_results/{c}/kma_index/masked_gene.comp.b",
        mb="output/4_cluster_results/{c}/kma_index/masked_gene.length.b",
        mc="output/4_cluster_results/{c}/kma_index/masked_gene.name",
        f="output/4_cluster_results/{c}/kma_index/flanks_with_gene.seq.b",
        m="output/4_cluster_results/{c}/kma_index/masked_gene.seq.b",
        j="output/4_cluster_results/{c}/kma_index/just_gene.seq.b"
    output:
        f="output/4_cluster_results/{c}/{c}.flanks_with_gene_dist",
        m="output/4_cluster_results/{c}/{c}.masked_gene_dist",
        j="output/4_cluster_results/{c}/{c}.just_gene_dist"
    params:
        f_db=lambda wildcards: "output/4_cluster_results/" +  wildcards.c + "/kma_index/flanks_with_gene",
        m_db=lambda wildcards: "output/4_cluster_results/" +  wildcards.c + "/kma_index/masked_gene",
        j_db=lambda wildcards: "output/4_cluster_results/" +  wildcards.c + "/kma_index/just_gene",
        dist_measure=config["distance_measure"]
     
    conda: "environment.yaml"
    shell:
        "kma dist -t_db {params.f_db} -o {output.f} -d {params.dist_measure};"
        "kma dist -t_db {params.m_db} -o {output.m} -d {params.dist_measure};"
        "kma dist -t_db {params.j_db} -o {output.j} -d {params.dist_measure}"




def dist_input(wildcards):
    checkpoint_output = checkpoints.split_cd_hit_results.get(**wildcards).output[0]
    return expand("output/4_cluster_results/{c}/{c}.masked_gene_dist",
           c=glob_wildcards(os.path.join(checkpoint_output, "{c}.tsv")).c)





rule finito:
    input:
        dist_input
    output:
        config="output/99_config.yaml",
        txt="output/5_plots/completed.txt"
    conda: "environment.yaml"
    shell:
        '''
        cp config.yaml {output.config};
        rm -r output/3_define_clusters/cd_hit_per_cluster;
        echo plots_completed > {output.txt};
        Rscript bin/plot_gene_clusters_from_flankophile.R
        '''





