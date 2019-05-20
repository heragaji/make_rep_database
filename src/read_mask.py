#!/usr/bin/env python

import csv
import os, sys
from Bio import SeqIO

def main():
    fasta_in = sys.argv[1]
    paf_in = sys.argv[2]
    fasta_out = sys.argv[3]
    d = {}
    with open(paf_in, 'r') as maskf:
        mask_reader = csv.reader(maskf, delimiter='\t')
        for row in mask_reader:
            if row[5] in d:
                d[row[5]].append((int(row[7]),int(row[8])))
            else:
                d[row[5]] = [(int(row[7]),int(row[8]))]
    with open(fasta_out,"w") as outf:
        for record in SeqIO.parse(fasta_in, 'fasta'):
            if record.id in d:
                for start_end in d[record.id]:
                    aligned = str(record.seq)[start_end[0]:start_end[1]].lower()
                    masked_seq = record.seq[:start_end[0]] + aligned + record.seq[start_end[1]+1:]
                record.seq = masked_seq
            SeqIO.write(record,outf,"fasta")


if __name__ == '__main__':
    main()