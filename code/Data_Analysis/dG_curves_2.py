import pandas as pd
import fit_utils as fit
import numpy as np
import melt_utils as mu

df = pd.read_csv('~/nn/code/Data_Analysis/testing_output.csv')

T = np. arange(0,101,5)

def df_to_dG(array):
    array = np.array(array[2:],dtype=float)
    try:
        curve = fit.fit_melt_curve_no_sigma(T,array)[0]
        dG = fit.get_dG(37,curve[0],curve[1])
    except:
        return np.nan
    return dG

def df_to_nupack_dG(array):
    seq = array[1]
    return mu.get_dG(seq)

def df_to_nupack_ddG(array):
    seq = array[1]
    return mu.get_ddG(seq)    

def df_to_nupack_dG_bp(array):
    seq = array[1]
    return mu.get_dG_bp_pfunc(seq)


#Formatting Data
df.ix[:,2:23] = df.ix[:,2:23].apply(pd.to_numeric)

df['fit_dG']     = df.apply(df_to_dG,axis=1,raw=True)
df['nupack_dG']  = df.apply(df_to_nupack_dG,axis=1)
df['nupack_ddG'] = df.apply(df_to_nupack_ddG,axis=1)
df['nupack_dG_bp'] = df.apply(df_to_nupack_dG_bp,axis=1)
df.to_csv('~/nn/code/Data_Analysis/testing_output_dG.csv')
