#!/usr/bin/env python
# coding: utf-8
from gensim.models.doc2vec import Doc2Vec, TaggedDocument
import sys
import pickle
import logging
import argparse
import gzip
import numpy as np 
from concurrent.futures import ProcessPoolExecutor


parser=argparse.ArgumentParser(description = 'train doc2vec model using variable number of output dimensions')
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
required.add_argument('--corpusFile', help='file', dest='corpus_file' )
required.add_argument('--targetTxLsf', action = 'store', dest = 'target_tx_lsf')
required.add_argument('--targetTxIDs', action  ='store', dest = 'target_txid_file')
required.add_argument('--edim', help='size of kmer',  type=int )
required.add_argument('--wc', help='size of kmer',  type=int )
required.add_argument('--dm', help='size of kmer',  type=int )
required.add_argument('--trainedModel', help='output file for lsf kmers', dest='out_model_file')
required.add_argument('--outMatrix', action = 'store', dest='out_matrix_file')

args=parser.parse_args()




logging.basicConfig(format='%(levelname)s : %(message)s', level=logging.INFO)
logging.root.level = logging.INFO


NCORES = 64

#PV-DM == skipgram == `dm=0`
print('training')
if args.load_model:
    model=Doc2Vec.load(args.out_model_file)
else:
    model = Doc2Vec(corpus_file=args.corpus_file, dm= args.dm, vector_size=,args.edim, min_count=args.wc, epochs=15, seed=42, workers=NCORES)
    model.delete_temporary_training_data(keep_doctags_vectors=False, keep_inference=True)
    model.save(args.out_model_file)
print('saved')


## paralell process without loading into memory


def wrap_infer(word_vec, model=model):
    infvec = model.infer_vector(word_vec[1:])
    str_vec = [str(i) for i in infvec]
    return [word_vec[0]] + str_vec


with open(args.target_tx_lsf) as targ_tx, \
        gzip.open(args.out_matrix_file, 'wt+') as out_mat, \
        open(args.target_txid_file) as txid_file, \
        ProcessPoolExecutor(max_workers=NCORES) as proc_exec:
    line = True
    while line:
        #first read in chunks to run on multiple processes
        curr_chunk = []
        chunk_size = NCORES*100
        #chunk_size=10
        i = 0
        while i < chunk_size and line:  # checks for eof from exhausted generator
            i += 1
            txid = txid_file.readline().strip('\n')
            line = targ_tx.readline()
            word_vec = [txid] + line.strip('\n').split(' ')
            if line != '':  # dont add exhausted generator
                curr_chunk.append(word_vec)
        ## now run mult
        if curr_chunk != []:
            res = proc_exec.map(wrap_infer, curr_chunk)
            res=list(res)
            for vec in res:
                out_mat.write(','.join(vec) + '\n')

    

        



