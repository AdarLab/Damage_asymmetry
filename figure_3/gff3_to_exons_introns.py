#!/usr/bin/env python
#####################################################################################
# This script gets gff file and generate bed files of exons and introns coordinates. 
#####################################################################################
import csv
import os
import sys
import gzip
from collections import defaultdict
import argparse

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

parser = argparse.ArgumentParser(description ="""This script parses input
                gff3 file and generates a table of introns""")
parser.add_argument('-i', '--infile', help="input file; default STDIN")
parser.add_argument('-o', '--outfiles', help="output files; default STDOUT", nargs="+")
parser.add_argument('-z', '--gzip', help="input is gzip compressed", \
                    action = 'store_true')
args = parser.parse_args()

if args.infile:
    if args.gzip:
        gff3_file = gzip.open(args.infile, 'rt')
    else:
        gff3_file = open(args.infile, 'r')
else:
    gff3_file = sys.stdin

if args.outfiles[0]:
    introns_file = open(args.outfiles[0], 'w')
else:
    introns_file = sys.stdout
if args.outfiles[1]:
    exons_file = open(args.outfiles[1], 'w')
else:
    exons_file = sys.stdout

def get_exons(gff3_file):
    tbl = csv.reader(gff3_file, delimiter = '\t')
    exons_dict = defaultdict(list)
    for line in tbl:
        if not line[0].startswith('#'):
            [
                chrom,
                feat_source,
                feat_type,
                start,
                stop,
                score,
                strand,
                phase,
                attribs
            ] = line
            if (feat_type == "exon"
                and 'transcript_id' in attribs):
                new_attribs = process_attribs(attribs)
                if 'gene_id' in new_attribs:
                    gene_id = new_attribs['gene_id']
                elif 'gene' in new_attribs:
                    gene_id = new_attribs['gene']
                else:
                    gene_id = 'UNKNOWN'
                tx = new_attribs['transcript_id']
                start, stop = int(start), int(stop)
                exons_dict[(chrom, feat_source, strand, gene_id.replace('"', ''), tx.replace('"', ''))].append((start, stop))
            elif(feat_type == "exon"
                 and 'Parent=transcript' in attribs):
                attribs = attribs.split(';')
                new_attribs = {}
                for attrib in attribs:
                    attrib = attrib.split('=')
                    new_attribs[attrib[0]] = attrib[1]
                gene_id = new_attribs['Parent'].strip('transcript:')
                tx = new_attribs['Parent'].strip('transcript:')
                start, stop = int(start), int(stop)
                exons_dict[(chrom, feat_source, strand, gene_id, tx)].append((start, stop))
    gff3_file.close()
    return exons_dict

def process_attribs(attribs):
    new_attribs = {}
    attribs = list(filter(None, attribs.split('; '))) ## removes empty strings, needed because some gff3 lines have ";;"
    for attrib in attribs:
        k, v = attrib.split(' ')
        if k == 'Dbxref':
            xrefs = v.split(',')
            for xref in xrefs:
                terms = xref.split(':')
                new_attribs[terms[-2]] = terms[-1]
        else:
            new_attribs[k] = v
    return new_attribs

def get_introns(exons_dict):
    introns_dict = defaultdict(set)
    for tx_info, exons in exons_dict.items():
        if len(exons) > 1:
            [
                chrom,
                feat_source,
                strand,
                gene_id,
                tx
            ] = tx_info
            exons = sorted(exons)
            i = 0
            while (i + 1) < len(exons):
                introns = (exons[i][1] + 1, exons[i+1][0] - 1)
                i = i + 1
                introns_dict[(chrom, feat_source, strand, gene_id, tx)].add(introns)
    for k in introns_dict.keys():
        print(k)
    return introns_dict
    
def tabulate_exons(exons_dict, exons_file):
    tbl = csv.writer(exons_file,
                     delimiter = '\t',
                     lineterminator = os.linesep)
    tbl.writerow(['#chrom',
                  'source',
                  'feature',
                  'exon_start',
                  'exon_end',
                  'score',
                  'strand',
                  'gene_id',
                  'tx_acc',
                  'exon_num',
                  'exon_ct'])
    for tx_info, exons in exons_dict.items():
        [
            chrom,
            feat_source,
            strand,
            gene_id,
            tx_acc
        ] = tx_info
        if strand == '+':
            exons = sorted(exons)
        elif strand == '-':
            exons = sorted(exons, reverse = True)
        num_exons = len(exons)
        for exon in exons:
            tbl.writerow([chrom,
                         feat_source,
                         "exon",
                         exon[0], #exon start
                         exon[1], #exon end
                         ".",
                         strand,
                         gene_id,
                         tx_acc,
                         exons.index(exon) + 1,
                         num_exons])
    exons_file.close()


def tabulate_introns(introns_dict, introns_file):
    tbl = csv.writer(introns_file,
                     delimiter = '\t',
                     lineterminator = os.linesep)
    tbl.writerow(['#chrom',
                  'source',
                  'feature',
                  'intron_start',
                  'intron_end',
                  'score',
                  'strand',
                  'gene_id',
                  'tx_acc',
                  'intron_num',
                  'intron_ct'])
    for tx_info, introns in introns_dict.items():
        [
            chrom,
            feat_source,
            strand,
            gene_id,
            tx_acc
        ] = tx_info
        if strand == '+':
            introns = sorted(introns)
        elif strand == '-':
            introns = sorted(introns, reverse = True)
        num_introns = len(introns)
        for intron in introns:
            tbl.writerow([chrom,
                         feat_source,
                         "intron",
                         intron[0], #intron start
                         intron[1], #intron end
                         ".",
                         strand,
                         gene_id,
                         tx_acc,
                         introns.index(intron) + 1,
                         num_introns])
    introns_file.close()

exons_dict = get_exons(gff3_file)                 
introns_dict = get_introns(exons_dict)
tabulate_exons(exons_dict, exons_file)
tabulate_introns(introns_dict, introns_file)
