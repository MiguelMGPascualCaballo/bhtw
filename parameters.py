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
## USED TO COMPUTE BOUNDS ##
############################
REL_TOL_1D = RBF(10)**(-4)
ABS_TOL_1D = RBF(10)**(-6)
MAX_ITERATIONS = 2**18
MAX_SUBINTERVALS = 2**15


############################
###### PRINTING INFOR ######
############################
VERBOSE = 1      # Quantity of printed information, from 0 to 2. (o prints nothing while 2 the most).
PRINT_DIGITS = 20
PRINT_EACH = 512


############################
# CONSTANTS, DO NOT CHANGE #
############################
SQRT_1_DIV_3 = RBF(3)**RBF('0.5') / RBF(3)
TWO = RBF(2)
ONE_DIV_2 = 1 / TWO
ONE_DIV_135 = 1/RBF(135)
ZERO = RBF(0)
PI=RBF(pi)


# In[ ]:




