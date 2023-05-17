#!/usr/bin/python

## This script accompanyies Monika Mioduchowska, Edyta Konecka, Bartlomiej Goldyn, Tom Pinceel, 
## Luc Brendonck, Dunja Lukic, Lukasz Kaczmarek, Tadeusz Namiotko, Katarzyna Zajac, Tadeusz Zajac, 
## Jan P. Jastrzebski, Krzysztof Bartoszek, Wolbachia enigma: Playing peekaboo with a master manipulator .

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

# In line 90 you need to provide the name of the fasta file that contains the input target sequences with which to compare the 
# collection of sequences that are present in the *fastaq files that are contained in the current directory.
# In line 77 you can change the file extenstion to something different than fastaq.
# It is assumed that each set of sequences in the fastaq files is split into two files.
# These files need to have names that follow alphabetically after each other, e.g. FileName_001.fastaq and FileName_002.fastq.


from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
from os import listdir
from glob import glob
import csv
import os

def sort_csv_by_score(input_csv,col_by_sort,b_reverse,b_to_int):
##https://stackoverflow.com/questions/47886931/sort-csv-by-column-name
    with open(input_csv, 'r') as f_input:
    csv_input = csv.DictReader(f_input,delimiter=';')
    if b_to_int:
            data = sorted(csv_input, key=lambda row: (int(row[col_by_sort]), row['Seq_ID']),reverse=b_reverse)
        else:
            data = sorted(csv_input, key=lambda row: (row[col_by_sort], row['Seq_ID']),reverse=b_reverse)

    output_csv = "sorted_"+input_csv

    with open(output_csv, 'wb') as f_output:    
    csv_output = csv.DictWriter(f_output, fieldnames=csv_input.fieldnames,delimiter=';', lineterminator=os.linesep)
    csv_output.writeheader()
    csv_output.writerows(data)
    
    f_input.close()
    f_output.close()

def hamming_distance(s1, s2):
## https://en.wikipedia.org/wiki/Hamming_distance#History_and_applications
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))
    
def my_dists_calc_all(fastq_records,seq_target,len_target,outfile_p_dist,outfile_align_score):
    for seqobj_unknown in fastq_records:
    seq_unknown = seqobj_unknown.seq
    len_unknown = len(seq_unknown)
    seq1 = seq_target
    seq2 = seq_unknown
    
    if (outfile_align_score!=""):
        align_score = pairwise2.align.localxx(seq_target, seq_unknown,score_only=True)
        outfile_align_score.write(seqobj_unknown.id+";"+str(align_score)+";"+"\n")
    
    if (outfile_p_dist!=""):
        alignments = pairwise2.align.localxx(seq_target, seq_unknown,score_only=False)	
        seq1 = alignments[0][0] 
        seq2 = alignments[0][1] 
        p_dist_value = ((1.0)*hamming_distance(seq1,seq2))/len(seq1)	
        outfile_p_dist.write(seqobj_unknown.id+";"+str(p_dist_value)+";"+"\n")

runnum = "01"  ## if one wants to rerun multiple times and have a unique identifier for each run
do_align_score = True ## should sequences be sorted accoring to their alignment score
do_p_dist = True ## should sequences be sorted accoring to their p-distance score

all_fastqfiles = glob('*fastq')
all_fastqfiles.sort()
input_amplified_sequence_fastaq_file_names=[None]*(len(all_fastqfiles)/2)
input_amplified_sequence_outnames=[None]*(len(all_fastqfiles)/2)

for j in range(0,(len(all_fastqfiles)/2)):
    input_amplified_sequence_fastaq_file_names[j]=[all_fastqfiles[2*j],all_fastqfiles[2*j+1]]
    if (all_fastqfiles[2*j].split('_')[0]==all_fastqfiles[2*j+1].split('_')[0]):
    input_amplified_sequence_outnames[j]=all_fastqfiles[2*j].split('_')[0]
    else:
    print "Problem with files "+str(all_fastqfiles[2*j])+" and "+str(all_fastqfiles[2*j+1])+" for i="+str(j)+".\n"


input_target_sequence_fasta_file_names=[["Wolbachia_OTU_univ_DMEL.fas"]]*len(input_amplified_sequence_fastaq_file_names)

for i in range(0,len(input_amplified_sequence_fastaq_file_names)):
    for target_sequence in input_target_sequence_fasta_file_names[i]:
    seqobj_targets = list(SeqIO.parse(target_sequence, "fasta"))
    for seqobj_target in seqobj_targets:
        seq_target = seqobj_target.seq
        name_target = seqobj_target.id 
        len_target = len(seq_target)
	
        outfile_align_score = ""
        outfile_p_dist = ""
        if (do_align_score):
	outfile_name_align_score = "dists_from_"+name_target+"_align_score_"+input_amplified_sequence_outnames[i]+"_"+runnum+"_localscore.csv"
            outfile_name_best_align_score = "bestdists_from_"+name_target+"_align_score_"+input_amplified_sequence_outnames[i]+"_"+runnum+"_localscore.csv"
	outfile_align_score = open(outfile_name_align_score,"w") 
	outfile_align_score.write("Seq_ID;align_score;"+"\n")

        if (do_p_dist):
	outfile_name_p_dist = "dists_from_"+name_target+"_p_dist_"+input_amplified_sequence_outnames[i]+"_"+runnum+"_localscore.csv"	    
	outfile_name_best_p_dist = "bestdists_from_"+name_target+"_p_dist_"+input_amplified_sequence_outnames[i]+"_"+runnum+"_localscore.csv"
	outfile_p_dist = open(outfile_name_p_dist,"w") 
	outfile_p_dist.write("Seq_ID;p_dist;"+"\n")

        for amplified_sequence_file in input_amplified_sequence_fastaq_file_names[i] :
	fastq_parser_amplified_sequences = SeqIO.parse(amplified_sequence_file, "fastq") 
	my_dists_calc_all(fastq_parser_amplified_sequences,seq_target,len_target,outfile_p_dist,outfile_align_score)
    
        if (do_align_score):
	outfile_align_score.close()
	sort_csv_by_score(outfile_name_align_score,"align_score",True,True)

        if (do_p_dist):
	outfile_p_dist.close()
	sort_csv_by_score(outfile_name_p_dist,"p_dist",False,False)

