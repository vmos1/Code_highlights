import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import glob 
import sys
import pickle 
import scipy.special as sc
# from scipy.special import logsumexp as sc.logsumexp

## Import fitting modules
import gvar as gv
import lsqfit

def f_plot_meff(t,y,max=200):
    '''
    Effective mass plot:
    t -> t values of correlator, y -> correlator gvar: val(err)
    max = max number of correlators to plot
    '''
    
    c=dict_global[str(s)]['c'] # Warning global dict_global and string s needed

    lst=[]

    n=(y[2:]+y[:-2])/(2.0*y[1:-1])
    m=np.arccosh(n) /c
#     print(m.shape)
    lst.append(np.nan)
    lst.extend(m)
    lst.append(np.nan)
    
    
    df1=pd.DataFrame(columns=['t','coeff','meff'])
    df1['t']=t
    df1['coeff']=y
    df1['meff']=np.array(lst)

    ## Plot effective mass
    plt.figure()
    df1=df1.dropna().head(max)
    val=df1.meff.values
    x,y,yerr=df1.t.values,[i.mean for i in val],[i.sdev for i in val]
    plt.errorbar(x,y,yerr,marker='o',linestyle='')


class corr:
    ''' Class to store correlators C(t) for specific l and s'''
    def __init__(self,df,l,c,s,Lt):
        self.df=df
        self.c=c
        self.s=s
        self.Lt=Lt
        self.l=l
        
        self.x_full=None
        self.y_full=None
        
        self.fit=None
        
    def f_perform_fit_exp(self,f_make_pars,func,fit_range,plt_range,num_exp=3,verbose=0,use_prior=False,plot=True):
              
        ## Define x and y for object for full plot range (Used only for fit)
        self.x_full={'t':self.df.t.values[min(plt_range):max(plt_range)+1],
                     'nterms':num_exp,'c':self.c,'Lt':self.Lt}
        self.y_full=self.df.coeff.values[min(plt_range):max(plt_range)+1]
        
        x=self.x_full.copy(); x['t']=x['t'][min(fit_range):max(fit_range)+1]
        y=self.y_full[min(fit_range):max(fit_range)+1]
        
        ### Create parameters dictionary: if prior gives gvars
        p0=f_make_pars(num_exp,self.l,use_prior)
        
        if use_prior:
            self.fit = lsqfit.nonlinear_fit(data=(x,y), prior=p0, fcn=func)
        else: 
            self.fit = lsqfit.nonlinear_fit(data=(x,y), p0=p0, fcn=func)

        ### Print fit details before plot
        if verbose==2: print(self.fit.format(maxline=True))
        elif verbose: print(self.fit.format(maxline=verbose))
        
        if plot: self.f_fit_plot(min(plt_range),max(plt_range),error_band=True,semilog=True)
        
        # Print description at the end
        for k in self.fit.p.keys():
            if not use_prior:
                print(k,'\tInit',self.fit.p0[k],'\t---Final',self.fit.p[k])
            else : 
                print(k,'\tInit',self.fit.prior[k],'\t---Final',self.fit.p[k])    

        print("chi-sqr",self.fit.chi2/self.fit.dof)

        
    def f_fit_plot(self,start=0,end=None,error_band=True,semilog=True):
        '''
        Function for plotting data with the fit lines and error bands.
        For correlators, using a semi-log plot.
        full_data=True, plots the entire data and the best-fit in the fit region
        '''
        
        ## Create y_predictions for fit function
        x=self.x_full; y=self.y_full
        curve_t={}
        curve_x=self.x_full.copy()
        curve_t['t']=np.linspace(min(x['t']),max(x['t']),500)
        curve_x['t']=curve_t['t']; ## Modify curve_x['ts'] to the range required
        
        curve_y=self.fit.fcn(curve_x,self.fit.p)
        obs_fit=gv.mean(curve_y)
        err_fit=gv.sdev(curve_y)
        sigma=2.0 # Width for band

        #### Plots
        fig=plt.figure()
          
        # Plot all data points
        x=self.x_full; y=self.y_full
        plt.errorbar(x['t'],y=gv.mean(y),yerr=gv.sdev(y),linestyle='None',color='black',marker='*',markersize=4)

        # Plot the best fit line
        plt.plot(curve_t['t'],obs_fit,color='blue')

        if error_band: # Plot an error band around the best fit line
            plt.fill_between(curve_t['t'],obs_fit-sigma*err_fit,obs_fit+sigma*err_fit,color='yellow')

        ## Plot data points used in fit with different color
        x=self.fit.x; y=self.fit.y
        plt.errorbar(x['t'],y=gv.mean(y),yerr=gv.sdev(y),linestyle='None',color='red',marker='H',markersize=5)    
        
        if semilog: plt.yscale('log')
        plt.ylabel('C(t)')
        plt.xlabel('t')
        plt.title("Plot")
        plt.show()
        
# Fit functions
def f_multi_exp(x,p):
    
    t=x['t']
    num_exp=x['nterms']
    c=x['c']
    Lt=x['Lt']
    
#     g=np.sqrt(8*np.pi/(3**(0.5)*Nsites))
    g=1.0
    
    ans=0.0
    ans=p['const']
    for i in range(num_exp):
        ans+= p["a{}".format(i)]*np.exp(-p["E{}".format(i)]* g * c *t)
    
    return ans


def f_multi_exp2(x,p):
    
    t=x['t']
    num_exp=x['nterms']
    c=x['c']
    Lt=x['Lt']
#     g=np.sqrt(8*np.pi/(3**(0.5)*Nsites))
    g=1
    
#     ans=0.0
    ans=p['const']
    for i in range(num_exp):
        ans+= p["a{}".format(i)]*np.exp(-p["E{}".format(i)]* g * c * t)
        ans+= p["a{}".format(i)]*np.exp(-p["E{}".format(i)]* g * c * (Lt-t))
    
    return ans

def f_cosh(x,p):

    t=x['t']
    num_exp=x['nterms']
    c=x['c']
    Lt=x['Lt']

#     g=np.sqrt(8*np.pi/(3**(0.5)*Nsites))
    g=1

    ans=0.0
    ans=p['const']
    for i in range(num_exp):
        ans+= p["a{}".format(i)]*np.cosh(p["E{}".format(i)]* g * c * (Lt/2-t))

    return ans


def f_make_pars(num_exp,l,use_prior):
    
    par_dict={}
    par_dict = gv.BufferDict()         # any dictionary works

    if l==0: 
        vals_a=[1,1,1,1]; vals_e=[1,3,5,7]
    elif l==2: 
        vals_a=[1,1,1,1]; vals_e=[2,3,5,7] # 2 is Wrap around term
    elif l==4: 
        vals_a=[1e-3,1,1,1]; vals_e=[3,5,7,9] # 3 is SB term
    elif l==6: 
        vals_a=[1e-5,1,1e-5,1,1]; vals_e=[1,7,3,9,11] # 1 is SB term
    else: 
        raise SystemError
    w_factor = 0.3 # fraction of uncertainty for prior
    
    keys_a=["a{0}".format(i) for i in range(num_exp)]
    keys_e=["E{0}".format(i) for i in range(num_exp)]
        
    for count,key in enumerate(keys_e):
        if use_prior: 
            par_dict[key]=gv.gvar(vals_e[count],vals_e[count]*w_factor)
        else : 
            par_dict[key]=vals_e[count]
    
    for count,key in enumerate(keys_a):
        if use_prior: 
            par_dict[key]=gv.gvar(vals_a[count],vals_a[count]*w_factor)
        else : 
            par_dict[key]=vals_a[count]
    
    par_dict['const']=1.0 # No constant term for 2pt function
    
    return par_dict

def f_scale_4pt(row,s,dict_global):
    a1=row['coeff']/(dict_global[str(s)]['2pt']**2) # divide by phi_2pt ^ 2 
    return a1

def f_get_data_df(fname,dict_global):

    arr=np.loadtxt(fname)

    s= int(fname.split('/')[-1].split('_')[2].split('k')[-1].split('t')[0])
    Lt=int(fname.split('/')[-1].split('_')[2].split('k')[-1].split('t')[-1])
    print(s,Lt)
    Nsites=10*(s**2)+2
    # c=dict_global[str(s)]['c']

    ### Gather data from file
    cols=['l','t','coeff']
    df_data=pd.DataFrame(arr,columns=cols)

    ### Process data
    # df_data.apply(lambda row: f_scale_4pt(row,s),axis=1)
    df_data['coeff']=df_data.apply(lambda row: f_scale_4pt(row,s,dict_global),axis=1)
    
    # df_data['coeff']=np.array([gv.gvar(i,i*0.01) for i in df_data.coeff.values])
    df_data['coeff']=np.array([gv.gvar(i,max(i*1e-4,1e-13)) for i in df_data.coeff.values]) # Ensure error doesn't go below machine precision
    

    return df_data,Lt


def f_sort_dict_energies(dict1):
    '''
    Sort dictionary by 'E' values. Use the order to sort 'a' values
    
    '''
    
    dict2=dict1.copy()
    # Extract keys starting with 'E'
    e_keys=[key for key in dict1.keys() if key[0]=='E']

    # Drop nan values 
    d2=dict([(key,dict1[key]) for key in e_keys if not np.isnan(gv.mean(dict1[key]))])
    
    # Sorting a dictionary by value
    d3={key:value for key, value in sorted(d2.items(), key=lambda x: x[1])}
    
    # Use sorted dictionary to rename dictionary entries. Eg: d_new['E2']=d_old['E1']
    for count,key in enumerate(d3.keys()):
        k_old=key[1:]
        ## Set new dict element = old dict 
        dict2['E{0}'.format(count)]=dict1['E{0}'.format(k_old)]
        dict2['a{0}'.format(count)]=dict1['a{0}'.format(k_old)]
    
    return dict2


def f_fit_all_models(df,dict_global,l,s,Lt,tmin_max=50,tmax=60,num_exp_min=1,num_exp_max=8):
    ''' Module to implement model averaging for different nexp and tmin values for fixed tmax'''
    
    df_fits=pd.DataFrame([])
    fit_lst=[]
    IC='AIC'
    
    count=0
#     full_params=["a{0}".format(i) for i in range(num_exp_max)]+["E{0}".format(i) for i in range(num_exp)]
    for num_exp in range(num_exp_min,num_exp_max+1):
        model_params=["a{0}".format(i) for i in range(num_exp)]+["E{0}".format(i) for i in range(num_exp)]+['const']
        keys=['name','tmin','tmax','chi2_dof','prob','Q','fit_wt','normed_prob'] + model_params
        dict1=dict.fromkeys(keys)
        dict1['num_exp']=num_exp
        
        for tmin in np.arange(0,tmin_max,1):
            if (df[(df.t>tmin) & (df.t<tmax)].shape[0]) >1: 
                dict1['name']='trange_{0}-{1}'.format(tmin,tmax)
                dict1['tmin'],dict1['tmax']=tmin,tmax
#                 ft=f_perform_fit_multi_exp(df,f_make_pars,f_multi_exp2,[tmin,tmax],[0,80],l,num_exp,0,use_prior=False,plot=False)
                Corr=corr(df=df,l=l,c=dict_global[str(s)]['c'],s=s,Lt=Lt)
                Corr.f_perform_fit_exp(f_make_pars,f_multi_exp2,np.arange(tmin,tmax),[0,tmax],num_exp,0,False,False) # l=0,1,2,3
                ft=Corr.fit
                for pars in model_params:  dict1[pars]=ft.p[pars]
        
                ## Sort dictionary by Energy values
#                 print(dict1)
                dict1=f_sort_dict_energies(dict1)
                
#                 dict1['fit_pass']=False if f_fit_pass(ft) else True
                
                fitIC, prob = get_raw_model_prob(ft, return_IC=True, IC=IC, ND=None, N_cut=tmin, yraw=None)
                for keys,vals in zip(['fit_wt','prob','Q','chi2_dof','dof','chi2','num_params','Ncut'],[fitIC,prob,ft.Q,(ft.chi2/ft.dof),ft.dof,ft.chi2,len(ft.p),tmin]):
                    dict1[keys]=vals
                
                df_fits=pd.concat([df_fits,pd.DataFrame(dict1,index=[count])])
                fit_lst.append(ft)
                count+=1
                
    return df_fits,fit_lst

############################
# Modules for model averaging 

def naive_IC_from_fit(fr):
    return fr.chi2

def BIC_from_fit(fr, ND):
    if ND is None:
        raise ValueError("Must specify ND to use Bayes IC!")
    return fr.chi2 + len(fr.p.keys()) * np.log(ND)

def AIC_from_fit(fr):
    return fr.chi2 + 2 * len(fr.p.keys())

def GAP_from_fit(fr):
    par_list = fr.p.keys()
    par_values = [fr.p[par] for par in par_list]
    Sigma_star = gv.evalcov(par_values)

    prior_values = [fr.prior[par] for par in par_list]
    Sigma_tilde = gv.evalcov(prior_values)

    return (
        fr.chi2
        - np.log(np.linalg.det(Sigma_star))
        + np.log(np.linalg.det(Sigma_tilde))
        + 2 * len(fr.p.keys())
    )

def get_raw_model_prob(fr, IC="AIC", N_cut=0, return_IC=False, ND=None, yraw=None):
    """
    Compute probability from log likelihood (LL) for a given fit.
    This is "raw" in the sense that it's not normalized over the
    space of all models, which should be done separately.

    Relation to info criteria:
        LL = -1/2 * IC

    Args:
      fr: Fit result object, from lsqfit module.
      IC: Which info criterion to use.  Options: AIC (default), BIC, GAP, BPIC, naive.
      N_cut: Correction for data subset selection.
      gamma_prior: Function with calculates the contribution to LL
                   from the presence of an indicator variable 'gamma'.
                   If 'None' (default), no gamma prior is included.
    """

    if IC == "BIC":
        LL = -0.5 * BIC_from_fit(fr, ND)
    elif IC == "AIC":
        LL = -0.5 * AIC_from_fit(fr)
    elif IC == "GAP":
        LL = -0.5 * GAP_from_fit(fr)
    elif IC == "naive":
        LL = -0.5 * naive_IC_from_fit(fr)
    else:
        raise ValueError(f"Unrecognized choice of info criterion: {IC}")

    # Correction to IC is +2 * N_cut - except for naive IC, which ignores this
    if IC != "naive":
        LL -= N_cut
        
    if return_IC:
        return LL, np.exp(LL)
    else:
        return np.exp(LL)


    
def model_avg(gv_list, pr_list):
    """
    Given a list of single-model expectation values {<f(a)>_M} as gvars,
    and a list of raw model probabilities, return the model-averaged estimate
    for <f(a)> as a gvar.
    """
#     print(gv_list)
    mean_avg = np.sum(gv.mean(gv_list) * pr_list)
    var_avg = np.sum(gv.var(gv_list) * pr_list)
    var_avg += np.sum(gv.mean(gv_list) ** 2 * pr_list)
    var_avg -= (np.sum(gv.mean(gv_list) * pr_list)) ** 2

#     print(np.sum(pr_list))
#     print(mean_avg,np.sqrt(var_avg))
#     return gv.gvar(mean_avg, np.sqrt(var_avg))
#     Divide by sum of probs, since sum don't add to 1 for subset
    return gv.gvar(mean_avg/np.sum(pr_list), np.sqrt(var_avg)/np.sum(pr_list)) 


def f_avg_all_models(df_fits):
    ''' Module to implement model averaging for different nexp and tmin values for fixed tmax'''
    
    ## Dataframe storing model averaged parameters and other info
    
    num_exp_max=np.max(np.unique(df_fits.num_exp.values))
    
    
    model_params=["a{0}".format(i) for i in range(num_exp_max)]+["E{0}".format(i) for i in range(num_exp_max)] + ['const']
    df_avg=pd.DataFrame(columns=model_params)
    
    # Cut data
    df_fits=df_fits[(df_fits.chi2_dof<2.0)]
#     df_fits=df_fits[df_fits.fit_pass]

    if df_fits.shape[0]<1: ## Fitting procedure didn't work well.
        for par in model_params:
            df_avg.loc['num_avg',par]=0
            df_avg.loc['value',par]=np.nan
        return df_fits,df_avg
    
    ### Normed prob for entire fit
    S=sc.logsumexp(df_fits.fit_wt.values)
#     print(S)
    df_fits['normed_prob']=df_fits.apply(lambda row : np.exp(row.fit_wt-S),axis=1)

    ## Select most probable fits and repeat normalization
    df_fits=df_fits[df_fits.normed_prob>1e-2]
    S=sc.logsumexp(df_fits.fit_wt.values)
    df_fits['normed_prob']=df_fits.apply(lambda row : np.exp(row.fit_wt-S),axis=1)

    
    # Compute model average parameters
    dict1={'name':'model_avg','fit_IC':None,'Q':None,'prob_AIC':None}
    
    for par in model_params:
        df_1=df_fits[['fit_wt',par]].dropna()
        # Normed prob for each parameter
        S=sc.logsumexp(df_1.fit_wt.values)
        df_1['normed_prob']=df_1.apply(lambda row : np.exp(row.fit_wt-S),axis=1)
#         print(par,df_1.shape)
        df_avg.loc['num_in_avg',par]=df_1.shape[0]
#         display(df_1)
        df_1=df_1[df_1.normed_prob>1e-2]
        
        dict1[par]=model_avg(df_1[par].values, df_1['normed_prob'].values)
        df_avg.loc['value',par]=dict1[par]

    dict1['normed_prob']=sum(df_fits.normed_prob)
#     print(dict1)
    df_fits=pd.concat([df_fits,pd.DataFrame(dict1,index=[df_fits.shape[0]+1])])

    return df_fits,df_avg
