
import numpy as np
import pandas as pd
import gvar as gv
import matplotlib.pyplot as plt
import itertools


def f_get_run_info_dict(dict1,input_dict,run_label):
    '''   
    Get dictionary with run info that will be added to Dataframe as columns
    '''
    keys=['Lx','Lt','beta','mf']
    for key in keys:
        dict1[key]=input_dict[key]
    
    run_key='beta-%s_mf-%s_Lx-%s_Lt-%s'%(input_dict['beta'],input_dict['mf'],input_dict['Lx'],input_dict['Lt'])
    
    if input_dict['F_action'] =='Mobius_dwf': ## Add Ls in label for DWF
        keys.append('dwf_Ls')
    
        
    run_key='beta-%s_mf-%s_Lx-%s_Lt-%s'%(input_dict['beta'],input_dict['mf'],input_dict['Lx'],input_dict['Lt'])
            
    if input_dict['F_action'] =='Mobius_dwf': ## Add Ls in label for DWF
        run_key+='_Ls-%s'%(input_dict['dwf_Ls'])

    dict1['run_label'] = run_label
    run_key+='_'+run_label 
    
    dict1['run_key']   = run_key
        
    return dict1

def f_conv_strg_to_gvar(x):
    '''
    convert entry from string to gvar
    '''
    try :
        ans=gv.gvar(x)
    except Exception as e:
        ans=gv.gvar(np.nan,np.nan)
        
    return ans


def f_plot_vary_beta(df,obs,ylabel='y',op_fname=0): 
    
    plt.figure()

    for Lx in np.unique(df['Lx'].values):
        for Lt in np.unique(df['Lt'].values):
            for mf in np.unique(df['mf'].values):
                for (run_label,marker) in zip(np.unique(df['run_label'].values), itertools.cycle('o*d>^sDHPdpx_')):
                    df1=df[(df.Lx==Lx)&(df.Lt==Lt)&(df.mf==mf)&(df.run_label==run_label)].sort_values(by=['beta'])
                    label=r'${%s}^3 \times %s ,\ m_f=%s,\ %s$'%(Lx,Lt,mf,run_label)

                    x=df1.beta.values
                    y=df1[obs].values
                    plt.errorbar(x,gv.mean(y),gv.sdev(y),linestyle='',label=label,marker=marker)



    beta_list=np.unique(df.beta.values)
    plt.show()
    plt.xlabel(r'$\beta$')
    plt.xticks(np.arange(min(beta_list),max(beta_list)+0.2,0.5))
    plt.legend()
    plt.ylabel(ylabel,rotation='vertical')
    # plt.xlim(10.2,11.2)
    # plt.ylim(0,0.04)

    if op_fname:
        plt.savefig(op_fname)
        
def f_plot_hmc_runs(df_input,column):
    ''' Plot behavior of quantity in MC time
    col = Plaquette, Polyakov, Traj_time, Accept'''
    
    keys_list=np.unique(df_input.run_key.values)
    assert len(keys_list)>0 ,"Shortened list has 0 elements" 
    
    plt.figure()

    
    for key,marker in zip(keys_list,itertools.cycle('>^*sDHPdpx_')):

        df=df_input[df_input.run_key==key]
        label=key

        x=df.iter.values
        if column=='Polyakov':
            y=np.abs(df[column].values)
        else:
            y=df[column].values

        plt.plot(x,y,linestyle='',label=label, marker=marker)

    plt.legend(loc='best')
    plt.ylabel(column)
    plt.xlabel('Trajectory')
    plt.show()

def f_scatter_plot(df_input,equil=700):
    '''
    Scatter plot of complex values of Polyakov loop with equilibriation
    '''
    
    keys_list=np.unique(df_input.run_key.values)
    assert len(keys_list)>0 ,"Shortened list has 0 elements" 
    
    fig=plt.figure()

    for key,marker in zip(keys_list,itertools.cycle('*s>^DHPdpx_')):
        
        df=df_input[df_input.run_key==key]
#         display(df)

        beta=np.unique(df.beta.values)[0]
        run_label=np.unique(df.run_label.values)[0]
        label='beta-%s_%s'%(beta,run_label)
        
        x=df.iter.values[equil:]
        y=df.Polyakov.values[equil:]

        if y.shape[0]<1: 
            print("%s points for %s"%(y.shape[0],label))
            
        y1=y.real
        y2=y.imag
        
        plt.scatter(y1,y2,label=label,marker=marker)

    plt.legend(loc='best')
    plt.xlabel('Real Polyakov loop')
    plt.ylabel('Imag Polyakov loop')
    plt.title("Scatter plot")
    
# Histogram

def f_ploop_histogram(df_input,equil=700,bins=100):
    '''
    Histogram of the magnitude of the Polyakov loop 
    '''

    keys_list=np.unique(df_input.run_key.values)
    assert len(keys_list)>0 ,"Shortened list has 0 elements" 
    
    fig=plt.figure()
    
    for key in keys_list:
        
        df=df_input[df_input.run_key==key]
#         display(df)

        beta=np.unique(df.beta.values)[0]
        run_label=np.unique(df.run_label.values)[0]
        label='beta-%s_%s'%(beta,run_label)

        x=df.iter.values[equil:]
        y=np.abs(df.Polyakov.values)[equil:]

        plt.hist(y,bins=100,label=label)
        plt.legend(loc='best')
#         plt.title("Polyakov loop Histogram for beta: %s"%(beta))


