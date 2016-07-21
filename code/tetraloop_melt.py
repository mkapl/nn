import melt_utils as mu
#import numpy as np
import itertools as it
stemLeft  = 'TCAGTGA'
stemRight = 'TCACTGA'

bases = ['A','C','G','U']
base = ['ACGU']
tetraloop = {}

i = 1
for p in it.product(bases,bases,bases,bases):
    tetraloop[''.join(p)] = [p]

for key in tetraloop:
    tetraloop[key] = [mu.get_melting_curves(stemLeft + key + stemRight)]
