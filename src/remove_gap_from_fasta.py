from Bio import SeqIO
import sys

fasta_in = sys.argv[1]
for record in SeqIO.parse(fasta_in, 'fasta'):
    print(">"+record.description)
    print(str(record.seq).replace('-',''))