import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ipywidgets import *
import glob 
import sys
import pickle 
import scipy.special as sc
# from scipy.special import logsumexp as sc.logsumexp

## Import fitting modules
import gvar as gv
import lsqfit

from modules import *


if __name__=="__main__":
    
    ### Set values for 2pt function and speed of light #####
    dict_global={}
    dict_global['4']=  {'2pt':1.126500594310e-02, 'c':2.758810226094e-01}
    dict_global['8']=  {'2pt':5.591911340388e-03, 'c':1.395714549742e-01}
    dict_global['16']= {'2pt':2.790831203460e-03, 'c':6.999296068404e-02}
    print(dict_global)
    
    ### Load entire data ####
    data_dir='../../data/free_theory/'

    data_dict={}
    s_list=[4,8,16]

    for s in s_list:
        print(data_dir+'s2xr_free_q5k{0}t*_4pt_pl.dat'.format(s))
        fname=glob.glob(data_dir+'s2xr_free_q5k{0}t*_4pt_pl.dat'.format(s))[0]
        print(fname,s)
        df,Lt=f_get_data_df(fname,dict_global)
        data_dict[str(s)]={'df':df,'Lt':Lt}
    
    ### Perform all fits and save results ###
    
    save_loc='../../data/stored_results/'

    cols=['s','l','E0','a0','E1','a1','const']
    df_results=pd.DataFrame(columns=cols)
    
    count=0 
    for s in [str(i) for i in s_list][:]:
        for l in [0,2,4,6][:2]:
            print(s,l)
            df_temp=data_dict[s]['df']
            df_temp=df_temp[df_temp.l==l][['t','coeff']]
            Lt=data_dict[s]['Lt']

            tmax=int(np.max(df_temp.t.values))
            print(df_temp.shape)        
            df_fits,fits_list=f_fit_all_models(df_temp,dict_global,l,s,Lt,tmin_max=int(Lt//4),tmax=tmax,num_exp_min=3,num_exp_max=4)

            ##### Save fits ########
            ## Save fit results in dataframe
            fname=save_loc+'l{0}_s{1}.df_fits'.format(l,s)
            print(fname)
            df_fits.to_pickle(fname)


    



