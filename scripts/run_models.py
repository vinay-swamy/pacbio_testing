#%%
#import tensorflow as tf 
import argparse 
import os 
import pathlib 
import numpy as np 
from loaders import *
from sklearn.ensemble import RandomForestClassifier 
from sklearn.model_selection import RandomizedSearchCV
from xgboost import XGBClassifier 
from functools import partial 
import pandas as pd 
from class_defs import Experiment
import pickle
#tf.test.is_gpu_available()
os.chdir('/data/swamyvs/pacbio_testing') 

#%%
def train_sk_model(model, dat):
    model.fit(dat.X_train, dat.Y_train)
    Y_pred_class = model.predict(dat.X_test)
    Y_pred_prob = model.predict_proba(dat.X_test)[:, 1]
    Y_true = dat.Y_test
    return model, Y_true,  Y_pred_class, Y_pred_prob


def parse_exp_sel(s):
    l=s.split(';')
    return [tuple(i.split(',')) for i in l]



# %%
default_lab_file = '/data/swamyvs/pacbio_testing/data/gtf_info/all_RPE_loose_target_tx.tsv'
parser=argparse.ArgumentParser()
parser.add_argument('--workingDir', action = 'store', dest = 'working_dir')
parser.add_argument('--inputFile', action = 'store', dest = 'input_file')
parser.add_argument('--labFile', action = 'store', dest='lab_file')
parser.add_argument('--inputType',action = 'store', dest='input_type')
parser.add_argument('--models', action = 'store', dest = 'models')
parser.add_argument('--expSelection', action='store', dest = 'exp_sel', default='all')
parser.add_argument('--nproc', action = 'store', type=int, default=2)
parser.add_argument('--saveDl', action = 'store_true', dest='save_dl')
parser.add_argument('--outputDir', action = 'store',  dest='output_dir')
args=parser.parse_args()
os.chdir(args.working_dir)
outdir= args.output_dir
lab_file = args.lab_file
if outdir[-1] is not '/':
    outdir+='/'
pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
#%%
####MODELS####
rf=partial(train_sk_model, RandomForestClassifier(n_estimators=100, random_state=32, n_jobs=args.nproc))
rf_weighted =partial(train_sk_model, RandomForestClassifier(n_estimators=100, random_state=32, n_jobs=args.nproc, class_weight='balanced'))

model_d={'rf' : rf, 'rf_weighted': rf_weighted}
##############




mode='file'

models=args.models.split(',')
if mode == 'file':
    with open(outdir+'model_results.csv', 'w+') as model_res_file:
        if args.input_type == 'skl':
            data_list=load_snakemake_data(args.input_file, lab_file, run_selection=args.exp_sel)
        elif args.input_type == 'ilearn':
            data_list = load_ilearn_data(args.input_file, lab_file)
        elif args.input_type == 'pyfeat':
            data_list = load_pyfeat_data(args.input_file, lab_file)
        elif args.input_type == 'all_ilearn':
            data_list = load_all_ilearn(lab_file)
        elif args.input_type == 'all_pyfeat':
            data_list = load_all_pyfeat(lab_file)
        elif args.input_type == 'merge':
            data_list = load_merge(lab_file, args.input_file, run_selection = args.exp_sel)
        elif args.input_type == 'load':
            with open('testing/data_list.pck', 'rb') as p:
                data_list=pickle.load(p)
        if args.save_dl:
            with open('testing/data_list.pck', 'wb+') as xf:
                pickle.dump(data_list, xf,  protocol=4)

        
        c_models= select_if_in(models, model_d.keys())
        assert len(c_models) > 0   
        c_model_dict={}
        for key  in c_models:
            c_model_dict[key] = model_d[key]

        trial=Experiment(data_list, c_model_dict, outdir, save_probs=True)
        if args.exp_sel == 'all':
            print('running all')
            trial_results=trial.run_all()
        else:
            print('running sel')
            trial_results=trial.run_selection(parse_exp_sel(args.exp_sel))
        [model_res_file.write(line) for line in trial_results ]
        
elif mode == 'folder':
    #TODO
    print('not_done_yet')






