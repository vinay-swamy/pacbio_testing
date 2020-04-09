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
from class_defs import Experiment, plot_metrics
from models import denseModel
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

def train_tf_model(model, outdir, obj):
    batch_size=int(obj.X_train.shape[0] / 4)
    nepochs=200
    history=model.fit(obj.X_train, obj.Y_train, epochs=nepochs, batch_size=batch_size, 
            validation_data=(obj.X_val, obj.Y_val), verbose=0)
    plot_metrics(history, outdir)
    pred=[p[1] for p in  model.predict(obj.X_test)]
    pred_class = np.int64([p > .5 for p in pred])
    true_class=[int(i[1]) for i in obj.Y_test]
    return model, true_class, pred_class, pred

def tune_hp_rf(outdir, obj):
    model=RandomForestClassifier(random_state=34)
    pars={'bootstrap': [True, False],
            'max_depth': [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, None],
            'max_features': ['auto', 'sqrt'],
            'min_samples_leaf': [1, 2, 4],
            'min_samples_split': [2, 5, 10],
            'n_estimators': [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]
            }
    #default validation is 5foldCV
    search=RandomizedSearchCV(estimator=model,param_distributions=pars, n_iter=100, scoring='f1', n_jobs=32, refit=True)
    search_res=search.fit(obj.X_train, obj.Y_train)
    best_model=search_res.best_estimator_
    Y_pred_class = best_model.predict(obj.X_test)
    Y_pred_prob = best_model.predict_proba(obj.X_test)[:, 1]
    Y_true = obj.Y_test
    return best_model, Y_true,  Y_pred_class, Y_pred_prob





def parse_exp_sel(s):
    l=s.split(';')
    return [tuple(i.split(',')) for i in l]



# %%
parser=argparse.ArgumentParser()
parser.add_argument('--workingDir', action = 'store', dest = 'working_dir')
parser.add_argument('--inputFile', action = 'store', dest = 'input_file')
parser.add_argument('--inputType',action = 'store', dest='input_type', default='skl')
parser.add_argument('--models', action = 'store', dest = 'models')
parser.add_argument('--expSelection', action='store', dest = 'exp_sel', default='all')
parser.add_argument('--nproc', action = 'store', type=int, default=2)
parser.add_argument('--saveDl', action = 'store_true', dest='save_dl')
parser.add_argument('--outputDir', action = 'store',  dest='output_dir')
args=parser.parse_args()
os.chdir(args.working_dir)
outdir= args.output_dir
if outdir[-1] is not '/':
    outdir+='/'
pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
#%%
####MODELS####
rf=partial(train_sk_model, RandomForestClassifier(n_estimators=100, random_state=32, n_jobs=args.nproc))
xgb=partial(train_sk_model,XGBClassifier(n_estimators=args.nproc, random_state=42)  )
dense=partial(train_tf_model, denseModel(), outdir)
hp_rf=partial(tune_hp_rf, outdir)
model_d={'random_forest' : rf , 'xgbc': xgb, 'dense': dense, 'hp_rf':hp_rf}
##############





mode='file'
lab_file = '/data/swamyvs/pacbio_testing/data/gtf_info/all_RPE_loose_target_tx.tsv'
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
        elif args.input_type == 'tf_merge':
            data_list = load_merge_tf(lab_file, args.input_file)
        elif args.input_type == 'skl-bal':
            data_list = load_snakemake_data(args.input_file, lab_file, run_selection=args.exp_sel, balance_data=True)
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






