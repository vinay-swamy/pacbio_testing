#%%
import pandas as pd 
import numpy as np
import pathlib
import glob
from class_defs import SklDataObj

def split_data_into_subsets(X_df, labs, X_df_val, val_labs, data_name):
    stringtie_cases=['all', 'stringtie', 'stringtie-pacbio','stringtie-scallop']
    scallop_cases=['all', 'scallop', 'scallop-pacbio','stringtie-scallop']
    labdf_all = labs.drop(['is_ism'], axis=1)
    labdf_stringtie = labdf_all[labdf_all.intersection_case.isin(stringtie_cases)].reset_index(drop=True)
    labdf_scallop = labdf_all[labdf_all.intersection_case.isin(scallop_cases)].reset_index(drop=True)
    labdf_no_ism=labs[labs.is_ism == 'not_ism'].drop(['is_ism'], axis=1).reset_index(drop=True)

    labdf_fam = (val_labs[val_labs['family'].isin(["krueppel_c2h2-type_zinc-finger_protein_family","g-protein_coupled_receptor_1_family"  ])] 
            .assign(target_label=lambda x: np.where(x['family'] == "g-protein_coupled_receptor_1_family", 1, 0 ) )
            .reset_index(drop=True))

    obj_list=[SklDataObj(data_name+ '_dataset=all',X_df, labdf_all, 'transcript', 'not_transcript'),
                SklDataObj(data_name+ '_dataset=stringtie', X_df, labdf_stringtie, 'transcript', 'not_transcript'),
                SklDataObj(data_name+ '_dataset=scallop', X_df, labdf_scallop, 'transcript', 'not_transcript'),
                SklDataObj(data_name+ '_dataset=not-ism', X_df, labdf_no_ism,  'transcript', 'not_transcript'),
                SklDataObj(data_name+ '_dataset=fam', X_df_val, labdf_fam,  'gprotein', 'zincfinger')
        ]
    return obj_list



def load_snakemake_data(X_file_name, lab_file,out_dir):
    fn=X_file_name.split('/')[-1]
    dm=int(fn.split('_')[2].split('-')[1])
    wc=int(fn.split('_')[3].split('-')[1])
    kmer_size=int(fn.split('_')[4].split('-')[1])
    edim=int(fn.split('_')[5].split('-')[1].split('.')[0] )
    #print([dm, wc, kmer_size, edim])
    fstem='/'.join(X_file_name.split('/')[:-1]) +'/'
    val_xfile=fn.split('_')
    val_xfile= fstem+ '_'.join(val_xfile[:2])+'_validation-tx_' + '_'.join(val_xfile[2:])
    print(val_xfile)
    val_labfile = 'data/gtf_info/all_RPE_loose_validation_tx.tsv'
    #lab_file='data/gtf_info/all_RPE_loose_target_tx.tsv'
    data_name=f'dm-{dm}_wc-{wc}_kmer-{str(kmer_size)}_edim-{str(edim)}'
    pathlib.Path(out_dir).mkdir(exist_ok=True,parents=True)
    all_positive_cases=['all', 'stringtie-pacbio', 'scallop-pacbio']

    X_df=pd.read_csv(X_file_name,names=['transcript_id']+ list(range(edim)))
    labs=(pd
        .read_csv(lab_file, sep= '\t')
        .assign(target_label=lambda x: np.where(x['intersection_case'].isin(all_positive_cases),1,0 )))
    X_df_val=pd.read_csv(val_xfile,names=['transcript_id']+ list(range(edim)))
    val_labs=pd.read_csv(val_labfile, sep='\t')
    skl_objs= split_data_into_subsets(X_df, labs, X_df_val, val_labs, data_name)
    return skl_objs


def load_ilearn_data(X_file_name, lab_file,out_dir):
    
    fn=X_file_name.split('/')[-1]
    #print([dm, wc, kmer_size, edim])
    fstem=X_file_name.split('/')[0] +'/'
    val_xfile=fstem + fn.split('_')[0]+ '_val_clean.tsv'
    val_labfile = '/data/swamyvs/pacbio_testing/data/gtf_info/all_RPE_loose_validation_tx.tsv'
    #lab_file='data/gtf_info/all_RPE_loose_target_tx.tsv'
    data_name=fn.split('.')[0]
    pathlib.Path(out_dir).mkdir(exist_ok=True,parents=True)
    all_positive_cases=['all', 'stringtie-pacbio', 'scallop-pacbio']

    col_size=pd.read_csv(X_file_name, sep='\t', nrows=20).shape[1] -2
    X_df=pd.read_csv(X_file_name, sep='\t', names=['transcript_id', 'dummy']+ list(range(col_size))).drop(['dummy'], axis=1)
    labs=(pd
        .read_csv(lab_file, sep= '\t')
        .assign(target_label=lambda x: np.where(x['intersection_case'].isin(all_positive_cases),1,0 )))
    X_df_val=pd.read_csv(val_xfile, sep='\t', names=['transcript_id', 'dummy']+ list(range(col_size))).drop(['dummy'], axis=1)
    val_labs=pd.read_csv(val_labfile, sep='\t')
    skl_objs= split_data_into_subsets(X_df, labs, X_df_val, val_labs, data_name)
    return skl_objs

def load_pyfeat_data(X_file_name, lab_file,out_dir):
    fn=X_file_name.split('/')[-1]
    fstem=X_file_name.split('/')[0] +'/'
    val_labfile = '/data/swamyvs/pacbio_testing/data/gtf_info/all_RPE_loose_validation_tx.tsv'
    data_name=fn.split('_')[0]
    pathlib.Path(out_dir).mkdir(exist_ok=True,parents=True)
    all_positive_cases=['all', 'stringtie-pacbio', 'scallop-pacbio']
    col_size=pd.read_csv(X_file_name, nrows=20).shape[1] -1
    X_df=pd.read_csv(X_file_name, names = list(range(col_size)) +['transcript_id'] ).iloc[:,::-1]
    labs=(pd
        .read_csv(lab_file, sep= '\t')
        .assign(target_label=lambda x: np.where(x['intersection_case'].isin(all_positive_cases),1,0 )))
    val_labs=pd.read_csv(val_labfile, sep='\t')
    X_df_val=X_df[X_df.transcript_id.isin(val_labs['transcript_id'])]
    skl_objs= split_data_into_subsets(X_df, labs, X_df_val, val_labs, data_name)
    return skl_objs



#def load_all_ilearn(lab_file, out_dir):
#%%

#all_dfs=[pd.read_csv(i, header=None, sep='\t') for i in files ]
def read_ilearn(file):
    try:
        col_size=pd.read_csv(file, sep='\t', nrows=20).shape[1] -2
        X_df=pd.read_csv(file, sep='\t', names=['transcript_id', 'dummy']+ list(range(col_size))).drop(['dummy'], axis=1)
        return X_df
    except:
        pass
def read_pyfeat(file):
    try:
        col_size=pd.read_csv(file, nrows=20).shape[1] -1
        X_df=pd.read_csv(file, names = list(range(col_size)) +['transcript_id'] ).iloc[:,::-1]
        return X_df
    except: 
        pass


def merge_all_csvs(files,  ftype):
    if ftype == 'ilearn':
        all_dfs=[read_ilearn(i) for i in files ]
    else: 
        all_dfs=[read_pyfeat(i) for i in files ]
    #%%
    all_dfs=[i for i in all_dfs if i is not None]
    final_df=all_dfs[0]
    for df in all_dfs[1:]:
        final_df= pd.merge(final_df, df, how= 'inner', left_on='transcript_id', right_on='transcript_id')
    return final_df
 
def load_all_ilearn(lab_file,out_dir):
    val_labfile = '/data/swamyvs/pacbio_testing/data/gtf_info/all_RPE_loose_validation_tx.tsv'
    #lab_file='data/gtf_info/all_RPE_loose_target_tx.tsv'
    data_name='all-ilearn'
    pathlib.Path(out_dir).mkdir(exist_ok=True,parents=True)
    all_positive_cases=['all', 'stringtie-pacbio', 'scallop-pacbio']
    train_files = glob.glob('/data/swamyvs/pacbio_testing/nonseq_feature_testing/ilearn_clean/*.tsv_clean.tsv')
    val_files = glob.glob('/data/swamyvs/pacbio_testing/nonseq_feature_testing/ilearn_clean/*.tsv_val_clean.tsv')
    X_df=merge_all_csvs(train_files, 'ilearn')
    X_df_val = merge_all_csvs(val_files, 'ilearn')
    labs=(pd
        .read_csv(lab_file, sep= '\t')
        .assign(target_label=lambda x: np.where(x['intersection_case'].isin(all_positive_cases),1,0 )))
    val_labs=pd.read_csv(val_labfile, sep='\t')
    skl_objs= split_data_into_subsets(X_df, labs, X_df_val, val_labs, data_name)
    return skl_objs

def load_all_pyfeat(lab_file,out_dir):
    val_labfile = '/data/swamyvs/pacbio_testing/data/gtf_info/all_RPE_loose_validation_tx.tsv'
    data_name='all-pyfeat'
    pathlib.Path(out_dir).mkdir(exist_ok=True,parents=True)
    all_positive_cases=['all', 'stringtie-pacbio', 'scallop-pacbio']
    all_files= glob.glob('/data/swamyvs/pacbio_testing/nonseq_feature_testing/pyfeat_out/*_features.csv')
    X_df=merge_all_csvs(all_files, 'pyfeat')
    labs=(pd
        .read_csv(lab_file, sep= '\t')
        .assign(target_label=lambda x: np.where(x['intersection_case'].isin(all_positive_cases),1,0 )))
    val_labs=pd.read_csv(val_labfile, sep='\t')
    X_df_val=X_df[X_df.transcript_id.isin(val_labs['transcript_id'])]
    skl_objs= split_data_into_subsets(X_df, labs, X_df_val, val_labs, data_name)
    return skl_objs




# %%
