
#%%

import argparse 
import os 
import pathlib 
import numpy as np 
from loaders import *
from sklearn.ensemble import RandomForestClassifier 
from xgboost import XGBClassifier 
from functools import partial 
import pandas as pd 
from class_defs import Experiment 
#os.chdir('/data/swamyvs/pacbio_testing') 

#%%
def train_sk_model(model, dat):
    model.fit(dat.X_train, dat.Y_train)
    Y_pred_class = model.predict(dat.X_test)
    Y_pred_prob = model.predict_proba(dat.X_test)[:, 1]
    return model, Y_pred_class, Y_pred_prob



# %%
parser=argparse.ArgumentParser()
parser.add_argument('--workingDir', action = 'store', dest = 'working_dir')
parser.add_argument('--inputFile', action = 'store', dest = 'input_file')
parser.add_argument('--inputType',action = 'store', dest='input_type', default='skl')
parser.add_argument('--nproc', action = 'store', type=int, default=2)
parser.add_argument('--outputDir', action = 'store',  dest='output_dir')
args=parser.parse_args()
os.chdir(args.working_dir)
#%%
outdir= args.output_dir
pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
mode='file'
lab_file = '/data/swamyvs/pacbio_testing/data/gtf_info/all_RPE_loose_target_tx.tsv'
if mode == 'file':
    with open(outdir+'model_results.csv', 'w+') as model_res_file:
        if args.input_type == 'skl':
            data_list=load_snakemake_data(args.input_file, lab_file, outdir)
        elif args.input_type == 'ilearn':
            data_list = load_ilearn_data(args.input_file, lab_file, outdir)
        elif args.input_type == 'pyfeat':
            data_list = load_pyfeat_data(args.input_file, lab_file, outdir)
        elif args.input_type == 'all_ilearn':
            data_list = load_all_ilearn(lab_file, outdir)
        elif args.input_type == 'all_pyfeat':
            data_list = load_all_pyfeat(lab_file, outdir)
        rf=partial(train_sk_model, RandomForestClassifier(n_estimators=100, random_state=32, n_jobs=args.nproc))
        xgb=partial(train_sk_model,XGBClassifier(n_estimators=args.nproc, random_state=42)  )
        model_d={'random_forest' : rf , 'xgbc': xgb}
        trial=Experiment(data_list, model_d, outdir)
        trial_results=trial.run_all()
        [model_res_file.write(line) for line in trial_results ]
        
elif mode == 'folder':
    print('')






