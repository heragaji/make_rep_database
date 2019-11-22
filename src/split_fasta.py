from Bio import SeqIO
import sys



def main():
    fasta_in = sys.argv[1]
    output_dir = sys.argv[2]
    top = sys.argv[3]
    for record in SeqIO.parse(fasta_in, 'fasta'):
        with open(output_dir+"/"+record.id+'/top'+top+"_"+record.id+".fa","w") as f:
            f.write(">"+record.id+'\n')
            f.write(str(record.seq)+'\n')

if __name__ == '__main__':
    main()