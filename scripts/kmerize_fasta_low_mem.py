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
parser.add_argument('--targetTx', action = 'store', dest = 'input_tx_file')
parser.add_argument('--outPfx', action='store', dest = 'out_prefix')
args=parser.parse_args()

full_lsf_file = args.out_prefix + '_full.lsf' 
full_txids_file = args.out_prefix + '_full_txids.txt'
targ_lsf_file = args.out_prefix + '_targ.lsf'
targ_txids_file = args.out_prefix + '_targ_txids.txt'

kmer_dist=1
kmer_size=args.kmerSize





with open(args.infasta) as infa_file,\
    open(args.input_tx_file) as inp_tx_file, \
    open(full_lsf_file, 'w+') as full_lsf, \
    open(targ_lsf_file, 'w+') as targ_lsf, \
    open(full_txids_file, 'w+') as full_tx, \
    open(targ_txids_file , 'w+') as targ_tx:
    selected_tx=set([tx.split('\t')[0] for tx in inp_tx_file])
    fasta = SeqIO.parse(infa_file, format='fasta')
    for record in fasta:
        header=str(record.id)
        kmers = kmerize_sequence(str(record.seq), kmer_size, kmer_dist, 1)
        line=' '.join(kmers[0]) + '\n'
        ids=header + '\n'
        full_lsf.write(line)
        full_tx.write(ids)
        if header in selected_tx:
            targ_lsf.write(line)
            targ_tx.write(ids)

