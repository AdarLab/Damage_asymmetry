"""This script gets an human sequence fasta files
and generates a counting table and a ratio table for the frequency
of different motifs"""

from Bio import SeqIO
import sys
import traceback
import warnings
import argparse
import csv
import random
import operator as op

NUC_NR = ['A','T']

MOT = ['A','T',2*'A',2*'T', 3*'A', 3*'T', 4*'A', 4 *'T',
       5*'A', 5*'T', 6*'A', 6*'T', 7*'A', 7*'T',
       8*'A', 8*'T', 9*'A', 9*'T', 10*'A', 10*'T']


def Counting(seq):
    """This function gets a sequence and returns the amount of motifs"""

    #Scan the sequence, looking for motifs

    counting = {k: 0 for k in MOT} # Initialize the counting dictionary.
    # Scan all the motifs and find them in the sequence
    for motif in MOT:
        if len(seq) > len(motif):  # Check if the sequence is longer than the motif itself.
            for i in range(len(seq)-len(motif)+1):
                if i == 0:  # In case the motif is in the beginning of the sequence
                    # print("start: " + seq[i:i+len(motif)] + " next nuc: " + seq[i+len(motif)])
                    if seq[i:i+len(motif)] == motif and seq[i+len(motif)] != motif[0]:  # Check if the next nucleotide is in not part of the motif.
                        counting[motif] += 1
                elif i == len(seq)-len(motif):  # In case the motif is in the end of the sequence
         
                    if seq[i:i+len(motif)] == motif and seq[i-1] != motif[0]: # Check if the previuos nucleotide is in not part of the motif.
                        counting[motif] += 1
                elif len(seq) > len(motif)+1: # In case the motif is in the middle of the sequence.
                    # Check if the motif is not part of another motif (e.g. TT is in TTT).

                    if seq[i:i+len(motif)] == motif and seq[i+len(motif)] != motif[0] and seq[i-1] != motif[0]:
                        counting[motif] += 1
    for nuc_nr in NUC_NR:
        counting[nuc_nr+"_NR"] = seq.count(nuc_nr)

    return counting

def parse(record):
    """ This function Uses Biopython's parse function to process individual FASTA records"""

    #Extract individual parts of the FASTA record

    identifier = record.id #The sequence's Id
    sequence = record.seq #The sequence itself
    sequence = sequence.upper() #Turns all the nucleotides to upper case

    return identifier, sequence

def append_counting(dict):
    """This function returns the counting table rows"""
    row_c = []
    # for nuc in NUC: #Scans all the elements and adds it to the table.
    #     row_c.append(dict[nuc])
    for mot in MOT:
        row_c.append(dict[mot])
    for nuc_nr in NUC_NR :
        row_c.append(dict[nuc_nr + "_NR"])
    #     #row.extend([dict["AA_NR"], dict["TT_NR"], dict["CC_NR"], dict["GG_NR"]])
    return row_c


def main():

    path_to_directory = sys.argv[1]
    path_to_file = sys.argv[2] #The file's path
    name = sys.argv[3] #The file's name
    name = name[:name.find('.')] #The organism's name
    print(name)
    counting_table = (open(path_to_directory + "/counting/" + name + "_counting.csv", "wb+")) #A file for the counting table
    counting_writer = csv.writer(counting_table)
    counting_header_keys = ["Id", "length"] + MOT + ["A_NR", "T_NR"]
    counting_writer.writerow(counting_header_keys) #Writes the header into file

    with open(path_to_file, mode='r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            Id,seq = parse(record)
            c = Counting(str(seq))
            counting_row = append_counting(c)
            counting_writer.writerow([Id, len(seq) - seq.count('N')] + counting_row) #Remove all unkown nucleotides from the lengrh.
    handle.close()
    counting_table.close()

main()
