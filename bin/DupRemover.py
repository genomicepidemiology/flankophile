#!/usr/bin/env python3
###################################
# File      - DupRemover
# Modified  - Sat Mar 27 14:03:43 CET 2021
# Sign      - Abhijeet
###################################
version = ': version (1.0.3)'
citation = '''\n\nCitation: Singh, Abhijeet. 2020. DupRemover: A Simple Program to Remove Duplicate Sequences from Multi-Fasta File 
GitHub: https://github.com/abhijeetsingh1704/DupRemover; DOI: 10.13140/RG.2.2.23842.86724.'''
###################################
import sys
import datetime
import subprocess
from collections import defaultdict
import argparse
#
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("Biopython missing, attempting to install...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython>=1.78"])
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
###################################
parser = argparse.ArgumentParser(prog='DupRemover.py', formatter_class=argparse.RawTextHelpFormatter, 
    description='Removes duplicate sequences in multifasta file, and append fasta header to unique sequence' + citation)
parser.add_argument("-i", "--input", dest='input', required=True, type=str, help='input fasta file')
parser.add_argument("-o", "--output", dest='output', required=False, type=str, help='output fasta file (default: Uniq_<input_fasta_file>)')
parser.add_argument("-v", "--verbose", dest='verbose', metavar='Y/y or N/n', default='Y', type=str, help='print progress to the terminal (default: verbose)')
parser.add_argument('-V', '--version', action='version', version='%(prog)s'+ str(version))
args = parser.parse_args()
###################################
# Input
input_file = args.input
input_obj = open(input_file, 'r')
###################################
# output
if args.output is not None:
    out_var = args.output
    output_obj = open(out_var, 'w')
else:
    out_var = "Uniq_" + input_file
    output_obj = open(out_var, 'w')
###################################
# verbosity
if args.verbose is not None:
    verbosity=args.verbose
else:
    verbosity=args.verbose(default)
verbosity=verbosity.upper() 
###################################
# Print info
date = datetime.datetime.now()
print("[Program]\t: DupRemover")
print("[Date]\t\t: "+ date.strftime("%Y-%m-%d %H:%M:%S"))
print("[Input file]\t: "+ input_obj.name)
print("[Output file]\t: "+ output_obj.name)
###################################
# reading and parsing fasta file
uniq_seqs = defaultdict(list)
for qry in SeqIO.parse(input_obj, 'fasta'):
    # making keys and list from sequence and ids + descriptions
    uniq_seqs[str(qry.seq)].append(qry.description)
# making unique sequences
final_seq = (SeqRecord(Seq(seqi), id="|".join(accn), name='', description='') for seqi, accn in uniq_seqs.items())
# write output file
output_num = SeqIO.write(final_seq, output_obj, 'fasta')
###################################
# verbosity for console
if verbosity == 'Y':
    print('-' * 25)
    for seqi, accn in uniq_seqs.items():
        print(' =|= '.join(map(str, accn)) + '\n' + seqi + '\n')
    print('-' * 25)
###################################
# close object
input_obj.close()
###################################
# sequence counts
input_num = 0
for input_seqs in open(input_file, 'r'):
    if input_seqs.startswith(">"):
        input_num += 1
# print stat   
print("[input seq]\t:", input_num)
print("[Output seq]\t:", output_num)
###################################
# Number of duplicate sequences
Duplicate_num = input_num - output_num
print("[Duplicates]\t:", Duplicate_num)
###################################
# close file
#input_file.close()
output_obj.close()
###################################
# END OF SCRIPT
