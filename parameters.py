#!/usr/bin/env python
# coding: utf-8

# In[1]:


from sage.all import *

############################
##### TYPES PRECISION ######
############################
prec = 200
RBF = RealBallField(prec)
CBF = ComplexBallField(prec)


############################
## BY DEFECT ##
############################
REL_TOL = RBF(10)**(-4)
ABS_TOL = RBF(10)**(-6)
MAX_ITERATIONS = 2**18
MAX_SUBINTERVALS = 2**15


############################
###### PRINTING INFOR ######
############################
VERBOSE = 0             # Quantity of printed information, from 0 to 1. (0 prints the essential, set to 1 to see progress and details).
VERBOSE_COUNTER_J = 40  # for treatment of J = (N, n) blocks 
VERBOSE_COUNTER_I = 512 # for integrals and image computators
VERBOSE_COUNTER_N = 100 # for polynomial_roots and ...
PRINT_DIGITS = 20
SHOW_ABCDE = True


############################
# CONSTANTS, DO NOT CHANGE #
############################
SQRT_1_DIV_3 = RBF(3)**RBF('0.5') / RBF(3)
TWO = RBF(2)
ONE_DIV_2 = 1 / TWO
ONE_DIV_135 = 1/RBF(135)
ZERO = RBF(0)
PI=RBF(pi)
MPI=-PI
TWOPI = TWO*PI
ONE = RBF(1)
SQRT2 = TWO**ONE_DIV_2

# ODE domain parameters
LEFT_DOM = MPI
LEN_DOM = TWOPI



# In[ ]:




