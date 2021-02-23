#!/usr/bin/env python
# Python 3.7.3 
######################################################################################################################
# This script gets an human sequence fasta files and generates a counting table for the frequency of different motifs.

# Parameters:
# 1) path_to_directory - The path to the directory for the tables.
# 2) path_to_file - The FASTA file itself.
######################################################################################################################

from Bio import SeqIO
import sys
import traceback
import warnings
import argparse
import csv
import random
import operator as op

NUC = ['A','T','C','G']
MOT = ['AA','AT','AC','AG',
       'TT','TA','TC','TG',
       'CC','CA','CT','CG',
       'GG','GA','GT','GC']

def Counting(seq):
    """This function gets a sequence and returns the amount of motifs"""

    counting = {k: 0 for k in MOT} # Intilaize a dictionary for the motif counts.
    #Scan the sequence, search for motifs and update the dictionary. 
    for i in range(len(seq)-1):
        if seq[i:i+2] in counting: # Dinucleotides.
            counting[seq[i:i+2]] += 1
    for nuc in NUC:
        counting[nuc] = seq.count(nuc) # Nucleotides.

    return counting

def parse(record):
    """ This function Uses Biopython's parse function to process individual FASTA records"""

    # Extract individual parts of the FASTA record.
    identifier = record.id # The sequence's Id.
    sequence = record.seq # The sequence itself.
    sequence = sequence.upper() # Convert the sequence to upper case.

    return identifier, sequence

def append_counting(dict):
    """This function returns the counting table rows"""
    row_c = []
    for nuc in NUC: # Scans all the elements and adds it to the table.
        row_c.append(dict[nuc])
    for mot in MOT:
        row_c.append(dict[mot])

    return row_c


def main():
    print sys.argv[1:]
    path_to_directory = sys.argv[1]
    path_to_file = sys.argv[2] 
    name = sys.argv[3] #The file's name
    name = name[:name.find('.')] # The organism's name

    counting_table = (open(path_to_directory + "/counting/" + name + "_counting.csv", "wb+")) # A file for the counting table

    counting_writer = csv.writer(counting_table)
    counting_header_keys = ["Id", "length"] + NUC + MOT # The header of the table
    counting_writer.writerow(counting_header_keys) # Write the header into file
    
    # Scan all the sequences in the given FASTA file and calculate the number of occurrences of the motifs.
    with open(path_to_file, mode='r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            Id,seq = parse(record)
            c = Counting(str(seq))
            counting_row = append_counting(c)
            counting_writer.writerow([Id, len(seq) - seq.count('N') ] + counting_row) # Remove all unkown nucleotides from the lengrh.

    handle.close()
    counting_table.close()

main()
