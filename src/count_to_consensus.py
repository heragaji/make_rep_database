import csv
import sys
import random
from argparse import ArgumentParser

random.seed(0)

def get_option():
    argparser = ArgumentParser()
    argparser.add_argument('-s', '--seq', type=str,
                           default="",
                           help='reference seq name')
    argparser.add_argument('-f', '--file', type=str,
                           default="",
                           help='input csv file name')
    argparser.add_argument('-g', '--gap', type=bool,
                           default=False,
                           help='alter gaps to -')
    return argparser.parse_args()

args = get_option()
seq_name = args.seq
filename = args.file
gap = ''
if args.gap:
  gap = '-'
print(">" + seq_name)
with open(filename, 'r') as f:
    reader = csv.reader(f,delimiter='\t')
    consensus = ""
    for row in reader:
        max_count = 0
        max_let = []
        if int(row[0]) > max_count:
            max_let = [gap]
            max_count = int(row[0])
        elif int(row[0]) == max_count:
            max_let += [gap]
        if int(row[1]) > max_count:
            max_let = ['A']
            max_count = int(row[1])
        elif int(row[1]) == max_count:
            max_let += ['A']
        if int(row[2]) > max_count:
            max_let = ['C']
            max_count = int(row[2])
        elif int(row[2]) == max_count:
            max_let += ['C']
        if int(row[3]) > max_count:
            max_let = ['G']
            max_count = int(row[3])
        elif int(row[4]) == max_count:
            max_let += ['G']
        if int(row[4]) > max_count:
            max_let = ['T']
            max_count = int(row[4])
        elif int(row[4]) == max_count:
            max_let += ['T']
        if int(row[5]) > max_count:
            max_let = ['N']
            max_count = int(row[5])
        elif int(row[5]) == max_count:
            max_let += ['N']
        consensus += random.choice(max_let)
    print(consensus)