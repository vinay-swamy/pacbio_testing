from Bio import SeqIO
import pickle
import sys 
import argparse


def kmerize_sequence(seq,k, n, label):
    '''
    k is the length of the kmer, and n is the distance between each kmer
    '''
    res = [seq[i:(i+k)] for i in range(0, len(seq)-k+1, n)]
    
    return((res, label))

parser=argparse.ArgumentParser(description = 'turn fasta into an lsf formatted kmer file')
parser.add_argument('--infasta', help='file')
parser.add_argument('--kmerSize', help='size of kmer',  type=int )
parser.add_argument('--outLsfFile', help='output file for lsf kmers', dest='out_lsf_file')
parser.add_argument('--txIDfile', help='output file for names of kmers', dest='out_tx_file')

args=parser.parse_args()
kmer_dist=1
kmer_size=args.kmerSize




with open(args.infasta) as infa_file, open(args.out_lsf_file, 'w+')  as ofl_lsf, open(args.out_tx_file, 'w+') as ofl_tx:
    fasta = SeqIO.parse(infa_file, format='fasta')
    for record in fasta:
        kmers = kmerize_sequence(str(record.seq), kmer_size, kmer_dist, 1)
        ofl_lsf.write(' '.join(kmers[0]) + '\n')
        ofl_tx.write(record.id + '\n')
