#%%
from functools import partial
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve, auc
from keras.utils import to_categorical
import pickle
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
plt.ioff() 

'''
The goal of these classes is to contain experiments, and allow to plug and play differnet models and data, with out rewriting a lot of code; right now I seem to be rewriting a bunch of stuff 

TODO:
move summary member function from skobj to dataobj, this will require a lot of code changing else where

''' 

def plot_metrics(history, outdir):
    metrics =  ['loss', 'auc', 'accuracy']
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for n, metric in enumerate(metrics):
        name = metric.replace("_"," ").capitalize()
        plt.subplot(2,3,n+1)
        plt.plot(history.epoch,  history.history[metric], color=colors[0], label='Train')
        plt.plot(history.epoch, history.history['val_'+metric],
                 color=colors[0], linestyle="--", label='Val')
        plt.xlabel('Epoch')
        plt.ylabel(name)
        if metric == 'loss':
            plt.ylim([0, plt.ylim()[1]])
        elif metric == 'auc':
            plt.ylim([0.7,1])
        else:
            plt.ylim([0,1])

        plt.legend()
    plt.suptitle('training metrics')
    plt.savefig(outdir + f'training_metrics.png')
    plt.close()


def model_results(obj,Y_true, Y_pred_class, Y_pred_prob, data_name, outdir):
    Y_true=np.asarray(Y_true)
    Y_pred_class = np.asarray(Y_pred_class)
    #Y_true = obj.Y_test

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    mpl.rcParams['figure.figsize'] = (12, 10)

    fpr, tpr, thresholds = roc_curve(Y_true, Y_pred_prob)
    AUC = roc_auc_score(Y_true, Y_pred_class)
    # Plot ROC curve
    plt.subplot(2, 1, 1)    
    plt.plot(fpr, tpr, label='ROC curve (area = %0.3f)' % AUC)
    plt.plot([0, 1], [0, 1], 'k--')  # random predictions curve
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate or (1 - Specifity)')
    plt.ylabel('True Positive Rate or (Sensitivity)')
    plt.title('Receiver Operating Characteristic')
    plt.legend(loc="lower right")

    pre, rec, thresholds = precision_recall_curve(Y_true, Y_pred_prob)
    AUC = auc(rec, pre)
    plt.subplot(2, 1, 2)
    plt.plot(rec, pre, label=' Prec/Rec (area = %0.2f)' % (AUC))
    plt.plot([1, 1], [1, 1], 'r--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall plot')
    plt.legend(loc="lower right")

    plt.suptitle(f'{data_name} ROC and PRC')
    #plt.show()
    plt.savefig(outdir + data_name + '_roc_prc.png')
    plt.close()
    #labs=[obj.zero_label, obj]
    cr_dict = classification_report(
        y_pred=Y_pred_class, y_true=Y_true, output_dict=True)
    tm = ['precision', 'recall', 'f1-score']
    zero_line = [str(cr_dict['0'][key]) for key in tm]
    one_line = [str(cr_dict['1'][key]) for key in tm]
    total_accuracy = str(cr_dict['accuracy'])
    one_accuracy = sum( Y_true[Y_true == 1] == Y_pred_class[Y_true == 1] ) / sum (Y_true == 1)
    zero_accuracy = sum( Y_true[Y_true == 0] == Y_pred_class[Y_true == 0] ) / sum (Y_true == 0)
    
    
    
    out_cs_line = ','.join([data_name] + zero_line +
                           one_line + [total_accuracy, str(one_accuracy), str(zero_accuracy)] ) + '\n'
    return(out_cs_line)





class DataObj:
    def __init__(self,name, x_train,y_train,x_test, y_test,  x_val=None,y_val=None):
        self.X_train=x_train
        self.Y_train=y_train
        self.X_test=x_test
        self.Y_test=y_test
        self.X_val=x_val
        self.Y_val=y_val
        self.name=name


class SklDataObj:
    def __init__(self,name, X_df, labs, one_label, zero_label):
        X_df_labeled=pd.merge(left=labs, how='inner', right=X_df, left_on='transcript_id', right_on='transcript_id')
        #assert X_df.shape[0]  == X_df_labeled.shape[0] # make sure number of rows match 
        X_data=np.asarray(X_df_labeled.iloc[:,3:])#drop the first 3 columns
        Y_vec=np.asarray(X_df_labeled['target_label'])
        self.Y_origin=X_df_labeled.iloc[:,:3]
        self.vec_y=Y_vec
        X_train, X_test, Y_train, Y_test =train_test_split(X_data,Y_vec,test_size=.2, 
                                                                                  random_state=42,stratify=Y_vec)
        DataObj.__init__(self,name, X_train, Y_train, X_test, Y_test)
        self.one_label=one_label
        self.zero_label=zero_label
    def summary(self):
        tr_len=len(self.X_train)
        ts_len=len(self.X_test)

        print(f'Training size: {tr_len}\ntesting size: {ts_len}')
        print(f'{self.one_label} count: {np.count_nonzero(self.vec_y == 1)}\n{self.zero_label} count : {np.count_nonzero(self.vec_y == 0)}')


class TfDataObj:
    def __init__(self, name ,X_df, labs, one_label, zero_label):
        assert labs.shape[1] == 3 #
        X_df_labeled=pd.merge(left=labs, right=X_df, left_on='transcript_id', right_on='transcript_id')
        X_data=np.asarray(X_df_labeled.iloc[:,3:])#drop the first 3 columns
        Y_vec=np.asarray(X_df_labeled['target_label'])
        X_train, X_val, Y_train_labs, Y_val_labs= train_test_split(X_data,labs,test_size=.2, random_state=42, stratify=Y_vec)
        X_train, X_test, Y_train_labs, Y_test_labs=train_test_split(X_train,Y_train_labs,test_size=.2, 
                                                                                  random_state=42,stratify=Y_train_labs['target_label'])
        Y_val=to_categorical(Y_val_labs['target_label'])
        Y_train=to_categorical(Y_train_labs['target_label'])
        Y_test=to_categorical(Y_test_labs['target_label'])
        DataObj.__init__(self, name, X_train,Y_train,X_test, Y_test,  X_val ,Y_val)
        self.one_label=one_label
        self.zero_label=zero_label

class Experiment:
    def __init__(self, obj_list,model_dict, outdir, save_probs=False, save_cm=False):
        obj_dict={}
        for obj in obj_list:
            obj_dict[obj.name]=obj
        self.obj_dict=obj_dict
        self.model_dict=model_dict
        self.trained_model_dict={}
        self.outdir=outdir
        self.save_probs = save_probs
        self.save_cm = save_cm
    

    def run_model(self, obj_name, model_name):
        obj=self.obj_dict[obj_name]
        model=self.model_dict[model_name]
        trained_model,y_true, y_pred_class, y_pred_prob=model(obj)
        
        if self.save_probs:
            pd.DataFrame.from_dict( {'y_true':y_true, 'y_pred_class':y_pred_class, 'y_pred_prob': y_pred_prob}).to_csv(self.outdir + 'y_probs.csv', index= False) 
        model_res_line = model_results(obj, y_true, y_pred_class, y_pred_prob, f'{obj_name}_{model_name}', self.outdir)
        return([model_res_line])
    
    def run_all(self):
        ''' run all models vs all data sets'''
        all_model_reslines=[]
        for obj_name in self.obj_dict.keys():
            for model_name in self.model_dict.keys():
                model_res=self.run_model(obj_name, model_name)
                all_model_reslines.append(model_res[0])
        return all_model_reslines


    def run_selection(self, pair_list):
        all_model_reslines=[]
        for pair in pair_list:
            model_name=pair[0]
            obj_name=pair[1]
            print(f'running {model_name} model and {obj_name} data set')
            model_res=self.run_model(obj_name, model_name)
            all_model_reslines.append(model_res[0])
        return all_model_reslines

            

