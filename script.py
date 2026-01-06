#!/usr/bin/env python
# coding: utf-8

# In[19]:


from sage.all import *
from parameters import *
from lemmas import *
import time


# In[ ]:


""" LOADING DATA """

with open('tw_data.txt', 'r') as f:
    lines = f.readlines()

data = [list(map(RBF, line.split())) for line in lines]

speed = data[0][0]
coeffs = data[0][1:]


# In[ ]:


#########################################
# lemma:twapprox_nonzero_coeffs
# Check that no coefficient is zero (necessary to do the argument)
for coef in coeffs:
    if coef.contains_zero():
        raise ValueError("Invalid coefficients: one or more coefficients contain zero.")


# In[21]:


""" APPROXIMATED TRAVELING WAVE RESIDUE """ 


# In[23]:


bound_norm_xi_L2="8.00338e-14"

#########################################
# lemma:residue_exis_L2
# 1.2147871400011354 seconds

# --- timing ---
t0 = time.perf_counter()
lemma_norm_xi_L2(speed, coeffs, bound_norm_xi_L2)
t1 = time.perf_counter()
# --- end timing ---
print(t1 - t0)


# In[24]:


""" RESULTS ABOUT THE POLYNOMIAL """


# In[26]:


radii_roots = "1.e-41"
la_real = "0.142128757204011"
la_imag = "-0.0707987316896479"


# In[27]:


# LOADS POLYNOMIAL ROOT APPROXIMATION (ONLY THE ONES INSIDE THE UNIT DISC)
roots = []
with open("polynomial_zeros.txt", "r") as f:
    for line in f:
        line = line.strip()
        rr, ii = line.split()
        root = CBF(rr,ii)  
        roots.append(root)


# In[28]:


radii_roots = RBF(radii_roots)
lamb = CBF(RBF(la_real), RBF(la_imag))


# In[29]:


#########################################
# lemma:enclose_poly_roots
# 56.015989612998965 seconds

# --- timing ---
t0 = time.perf_counter()
statement_lemma_enclose_roots = lemma_enclose_roots(speed, coeffs, roots, radii_roots)
t1 = time.perf_counter()
# --- end timing ---

print(t1 - t0)


# In[35]:


#########################################
# lemma:Proots1
# 16.029408063999654 seconds

# --- timing ---
t0 = time.perf_counter()
statement_lemma_roots_compatible_exis = lemma_roots_compatible_exis(speed, coeffs, roots, radii_roots)
t1 = time.perf_counter()
# --- end timing ---

print(t1 - t0)


# In[31]:


#########################################
# lemma:Proots2
# 15.709324043000379 seconds

# --- timing ---
t0 = time.perf_counter()
statement_lemma_roots_compatible_stab = lemma_roots_compatible_stab(speed, coeffs, roots, radii_roots, lamb)
t1 = time.perf_counter()
# --- end timing ---

print(t1 - t0)


# In[32]:


""" NORMS FOR LINEAR INVERSE OPERATOR, EXISTENCE """


# In[33]:


bound_beta_maximum = "1.99056e1"

#########################################
# lemma:beta_max
# 59 iterations
# 1.3269569469994167 seconds

# --- timing ---
t0 = time.perf_counter()
lemma_beta_maximum(speed, coeffs, bound_beta_maximum)
t1 = time.perf_counter()
# --- end timing ---
print(t1 - t0)


# In[34]:


bound_betamod_maximum = "5.07194"

#########################################
# lemma:beta_mod_max_exis
# 103 iterations
# 15.027077650000138 seconds

# --- timing ---
t0 = time.perf_counter()
lemma_betamod_maximum(speed, coeffs, bound_betamod_maximum)
t1 = time.perf_counter()
# --- end timing ---
print(t1 - t0)


# In[37]:


bound_betamod2_L2sq = "1.21079e4"
abs_tol = '1e-2' 
rel_tol = '1e-4'

#########################################
# lemma:betamod2_exis_L2
# 3597 iterations
# 22 minutes and 38.70595273900108 secs seconds
# 0.01 as abs_tol and 0.0001 as rel_tol

# --- timing ---
t0 = time.perf_counter()
lemma_betamod2_L2sq(speed, coeffs, bound_betamod2_L2sq, abs_tol, rel_tol)
t1 = time.perf_counter()
# --- end timing ---
print(t1 - t0)


# In[36]:


bound_betaiota = "7.74418e4"
abs_tol = "1e0"

#########################################
# lemma:betaiota_primi_exis_L1
# 2577 iterations
# 1 hour 29 minutes and 5.604465806000007 seconds
# 1 as abs_tol

# --- timing ---
t0 = time.perf_counter()
lemma_betaiota(speed, coeffs, bound_betaiota, abs_tol)
t1 = time.perf_counter()
# --- end timing ---
print(t1 - t0)

