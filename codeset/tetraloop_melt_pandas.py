import melt_utils as mu
import itertools as it
import fit_utils as fit
import numpy as np
import pandas as pd



weak = ['A','U']
strong = ['G','C']
bases = ['A','C','G','U']
base = ['ACGU']
tetraloop = []
Trange = range(0, 101, 5)
nopen = 1 
CtoK = 273.15

def ReverseComplement1(seq):
    seq_dict = {'A':'U','U':'A','G':'C','C':'G'}
    return "".join([seq_dict[base] for base in reversed(seq)])

def get_Fluor(array,T,nopen):
    T = T + CtoK
    seq = array[0]
    prop = mu.get_prop_closing_pair(seq,T,nopen)
    return prop

for stem in it.product(strong,weak):
    for tetra in it.product(bases,bases):
        key = ''.join(stem)  + ''.join(tetra) + ReverseComplement1(''.join(stem))
        tetraloop.append(key)

series = pd.Series(tetraloop)
tetraloop_df = pd.Series.to_frame(series)

for T in Trange:
   col_name  = str(T) +'_' + str(nopen)
   tetraloop_df[col_name] = tetraloop_df.apply(get_Fluor,axis=1,args=(T,nopen))



tetraloop_df.to_csv('~/nn/code/Data_Analysis/SWWWW.csv')
#print tetraloop_df
#with open('full_4base_stem_WWWW_dG','w') as f:
#    f.write('%s,' % range(0,101,5))
#    f.write('dG')
#    f.write('\n')
#    for key in tetraloop:
#        f.write('%s,' % key)
#       # f.writelines(['%s' % item for item in tetraloop[key]])
#        f.writelines(','.join(map(str,tetraloop[key])))
#        f.writelines(','.join(map(str,mu.get_dG(tetraloop[key]))))
#        f.write('\n')
