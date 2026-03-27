#!/usr/bin/env python
# coding: utf-8

# In[7]:


from sage.all import *
from lemmas import residual_norms, non_degen_coeffs, integral_bounds, hk_norms, svd_results, ode_resis_results, substituting_estimates
from time import perf_counter
from printing_macros import print_lemma, print_time


# In[8]:


def run_block(block):
    for lemma in block:
        start = perf_counter()
        verified, lemma_label, bounds = lemma()
        elapsed = perf_counter() - start
        print_lemma(verified, lemma_label, bounds, elapsed)

def run_substitutions():
    verified, description = substituting_estimates()
    if verified:
        print(f"Every substitution is OK.")
    else:
        print(description)
        


# In[ ]:


blocks = [
    residual_norms,   # checked
    non_degen_coeffs, # checked
    integral_bounds,  # checked
    hk_norms,         # checked
    svd_results,      # checked 
    ode_resis_results # checked
]

for x in blocks:
    run_block(x)

run_substitutions()





