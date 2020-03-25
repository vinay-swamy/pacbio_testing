
#%%

import argparse 
import os 
import pathlib 
import numpy as np 
from loaders import load_ilearn_data, load_snakemake_data 
from sklearn.ensemble import RandomForestClassifier 
from xgboost import XGBClassifier 
from functools import partial 
import pandas as pd 
from class_defs import Experiment 
os.chdir('/data/swamyvs/pacbio_testing') 

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
#%%
outdir= args.output_dir
pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
mode='file'
if mode == 'file':
    with open(outdir+'model_results.csv', 'w+') as model_res_file:
        if args.input_type == 'skl':
            data_list=load_snakemake_data(args.input_file, 
                            'data/gtf_info/all_RPE_loose_target_tx.tsv',
                            outdir)
        elif args.input_type == 'ilearn':
            data_list = load_ilearn_data(args.input_file, 
                            'data/gtf_info/all_RPE_loose_target_tx.tsv',
                            outdir)
        rf=partial(train_sk_model, RandomForestClassifier(n_estimators=100, random_state=32, n_jobs=args.nproc))
        xgb=partial(train_sk_model,XGBClassifier(n_estimators=args.nproc, random_state=42)  )
        model_d={'random_forest' : rf , 'xgbc': xgb}
        trial=Experiment(data_list, model_d, outdir)
        trial_results=trial.run_all()
        [model_res_file.write(line) for line in trial_results ]
        
elif mode == 'folder':
    print('')






