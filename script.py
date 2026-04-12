#!/usr/bin/env python
# coding: utf-8

# In[1]:


from sage.all import *

from time import perf_counter
from printing_macros import print_lemma, print_subst, print_time
from lemmas import (
    residual_norms, 
    non_degen_coeffs, 
    integral_bounds, 
    hk_norms, 
    svd_results,
    ode_resis_results,
    substituting_estimates
)


# In[2]:


def run_block(block):
    for lemma in block:
        start = perf_counter()
        verified, lemma_label, bounds = lemma()
        elapsed = perf_counter() - start
        print_lemma(verified, lemma_label, bounds, elapsed)

def run_substitutions():
    verified, description = substituting_estimates()
    print_subst(verified, description)
        


# In[3]:


blocks = [
    residual_norms, 
    non_degen_coeffs,
    integral_bounds, 
    hk_norms,      
    svd_results, 
    ode_resis_results,
]

for x in blocks:
    run_block(x)
    
run_substitutions()


# In[ ]:




