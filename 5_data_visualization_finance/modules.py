import numpy as np
import pandas as pd

def f_interest(x0,r,n,mode='compound'): # interest
    '''
    Computes interest. 2 modes: simple and compound
    Args: 
        x0-> Initial amount
        r -> rate of interest in percent. eg: 0.01
        n -> Number of terms
        mode-> 'compound'(default) or 'simple'.
    Returns:
        Final amount
    '''
    df = pd.DataFrame(columns=['term','amount','gain'])
    arr = np.zeros((n+1,3))
    arr[0] = np.array([0,x0,0])

    if mode=='compound':
        amount=x0*((1+r)**n)

        # arr = np.array([(i,x0+i*r*x0,i*r*x0)for i in range(n+1)])
        
        for i in range(1,n+1):
            arr[i,0] = i
            arr[i,1] = arr[i-1,1] * (1+r) 
            arr[i,2] = arr[i,1]- x0
            
    elif mode=='simple':
        amount=x0+x0*r*n
        # arr = np.array([(i,x0+r*x0,r*x0)for i in range(n+1)])
        for i in range(1,n+1):
            arr[i,0] = i
            arr[i,1] = arr[i-1,1] + r * x0
            arr[i,2] = arr[i,1]- x0

    for idx,col in enumerate(df.columns):
        df[col] = arr[:,idx]

    return amount,df


# SIP using formula
def f_SIP(x0,r,n):
    '''
    Computes the SIP using formula. 
    ie. You invest a fixed amount x0 every term at an interest 'r' for 'n' terms.
    
    Arguments : 
            x0 -> principal,
            r -> rate of interest, 
            n -> number of terms.
            check_mode: True-> computes by brute force to check with formula. Default is False.
    '''       

    # Computing using formula 
    R=1+r
    if n==0:   amount=x0
    else:      amount=(x0*(R**(n+1)-1)/(R-1))
    
    ## Compute table
    df = pd.DataFrame(columns=['term','amount','gain'])
    arr = np.zeros((n+1,3))
    arr[0] = np.array([0,x0,0])
    
    for i in range(1,n+1):
        arr[i,0] = i
        arr[i,1] = arr[i-1,1] * (1+r) + x0 
        arr[i,2] = arr[i,1]- x0*(i+1)
    
    for idx,col in enumerate(df.columns):
        df[col] = arr[:,idx]
        
    return amount,df


def f_EMI(x0,r,X_or_n,mode='amount'):
    '''
    Computes the EMI using formula. 
    ie. You take a loan x0. You pay a fixed amount X every term at an interest 'r' for 'n' terms, so that 
    the amount is paid.
    
    One of 'X' and 'n' are known. The other can be found
    Arguments : 
            x0 -> principal,
            r -> rate of interest, 
            n -> number of terms.
            mode: 'amount' (computes X) or 'period' (computes n)
            check_mode: True-> computes by brute force to check with formula. Default is False.
    
    For 'amount': x0,r,n
    For 'period': x0,r,X
    
    Returns : X or n

    '''
    
    df = pd.DataFrame(columns=['term','amount','term payment','total payment'])
    
    # Computing using formula    
    R=1+r

    if mode=='amount':
        n=X_or_n # Period has been given
        assert n>0,"Number of terms has to be greater than zero"
        
        X=(x0*(R-1)/(1-(1.0/R**(n))))
        
        
        arr = np.zeros((n+1,4))
        arr[0] = np.array([0,x0,0,0])
        
        for i in range(1,n+1):

            arr[i,0] = i
            arr[i,1] = np.around(arr[i-1,1] * R - X,5 )
            arr[i,2] = np.around(X,5)
            arr[i,3] = np.around(X * i,5)
            if ((np.around(arr[i,1],2)==0) and i==n ):
                arr[i,1]=0.0

        for idx,col in enumerate(df.columns):
            df[col] = arr[:,idx]

        return X,df
    
    elif mode=='period':
        X=X_or_n # Amount has been given
        if (X < x0 * r) : 

            print("EMI {} is less than interest accrued {}. Not possible to repay. Need to use a higer EMI".format(X,x0*r))
            n = np.inf
        
        elif X >= x0:
            n=1 # Exception required. Formula will give log of -ve number.
        else:
            Nr=1.0/(1-(x0/X)*(R-1)) # Taking the -1 factor inside to avoid issue with logarithms.
            Nr=(1-(x0/np.float64(X))*(R-1)) # Taking the -1 factor inside to avoid issue with logarithms.
            n=-1.0*np.log(Nr)/np.log(R)

            n_int=int(np.ceil(n))
            arr = np.zeros((n_int+1,4))
            arr[0] = np.array([0,x0,0,0])

            for i in range(1,n_int+1):

                arr[i,0] = i
                if i!=n_int: 
                    arr[i,1] = np.around(arr[i-1,1] * R - X,5 )
                    arr[i,2] = np.around(X,5)
                    arr[i,3] = np.around(X * i,5)     
                elif (i==n_int):
                    arr[i,1] = 0.0
                    arr[i,2] = np.around(arr[i-1,1]*R ,5)   ### Last term, just pay the balance 
                    arr[i,3] = arr[i-1,3] + arr[i,2]         

            for idx,col in enumerate(df.columns):
                df[col] = arr[:,idx]

        return n,df
