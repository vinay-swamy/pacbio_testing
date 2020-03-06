#!/usr/bin/env python
# coding: utf-8


from gensim.models.doc2vec import Doc2Vec, TaggedDocument
import sys
import pickle
import logging
import argparse


parser=argparse.ArgumentParser(description = ' train doc2vec model using variable number of output dimensions')
parser.add_argument('--corpusFile', help='file', dest='corpus_file' )
parser.add_argument('--dimSize', help='size of kmer',  type=int )
parser.add_argument('--trainedModel', help='output file for lsf kmers', dest='outfile')

args=parser.parse_args()

ndims=args.dimSize


logging.basicConfig(format='%(levelname)s : %(message)s', level=logging.INFO)
logging.root.level = logging.INFO



#PV-DBOW == skipgram == `dm=0`
print('training')
model = Doc2Vec(corpus_file=args.corpus_file, dm=0, vector_size=ndims, min_count=3, epochs=15, seed=42, workers=64)





model.save(args.outfile)

