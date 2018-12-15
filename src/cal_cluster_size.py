import sys
from Bio import SeqIO

fasta_in = sys.argv[1]
for record in SeqIO.parse(fasta_in, 'fasta'):
    score = 0.01
    print(str(record.id)+"," + str(score))
