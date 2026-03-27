#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from sage.all import *
from sage.rings.real_arb import RealBall
from sage.modules.free_module_element import vector

from collections import deque

from parameters import RBF, CBF, VERBOSE
from printing_macros import print_RBF
from load_data import interval_from_J


# In[ ]:


def o(out, verbose=VERBOSE):
    """Verifies that out is RBF"""
    
    #######################################
    # Verifies that out is a RealBall (RBF)
    if not isinstance(out, RealBall):
        raise TypeError(f"out must be of type RealBall (RBF) but it is {type(out)}")

    #######################################
    # Verifies that no precision was lost during the computation
    # (i.e., the parent ring/type is still RBF)
    if out.parent() != RBF:
        raise TypeError(f"Precision was lost during the process, out is: {out.parent()}")

    if verbose:
        check_str(f'{print_RBF(out)}') 
    return True


# In[ ]:


def oo(v, n):
    if len(v) != n:
        raise ValueError(
            f"Length mismatch: v has {len(v)} entries, should be {n}"
        )
    
    for out in v:
        assert o(out, verbose=0)

    check_str(f'iterable of {n} RBF elements')
    return True


# In[ ]:


def square_mat_order(M, ord_mat, name_mat='M'):
    
    if M.nrows() != ord_mat or M.ncols() != ord_mat:
        raise ValueError(
            f"{name_mat} must be {ord_mat}x{ord_mat} but is {M.nrows()}x{M.ncols()}."
        )

    check_str(f'{name_mat} is a square matrix of order {ord_mat}')
    return True
    


# In[ ]:


def ball_vector(v):
    if not hasattr(v, "base_ring"):
        raise TypeError("Input must be a Sage vector")

    base = v.base_ring()

    if base is not RBF and base is not CBF:
        raise TypeError(
            f"Vector base ring must be RBF or CBF, got {base}"
        )
        

    check_str(f'vector with base {base}') 
    return True


# In[ ]:


def data_block(values, n):
    # Infer polynomial order from the first block
    J0, values0 = values[0]
    if len(values0) != n:
        raise ValueError(
            f"First values entry has len(values)={len(values0)} but expected {n} (=len(vs))."
        )

    poly_order = len(values0[0]) - 1
    if poly_order < 0:
        raise ValueError("Invalid polunomial data: values0[0] is empty, cannot infer poly_order.")

    #######################################
    # Structural checks for every (J, values)
    for counter, blk in enumerate(values):

        ###################################
        # Each blk must be a pair (J, values)
        if not isinstance(blk, (list, tuple)) or len(blk) != 2:
            raise ValueError(f"Malformed values blk at index {counter}: expected (J, values).")

        J, values = blk

        ###################################
        # values must have length n
        if len(values) != n:
            raise ValueError(
                f"Malformed values entry at index {counter}: "
                f"len(values) = {len(values)} but expected {n} (=len(vs))"
            )
    
        ###################################
        # Each values[idx] must be a list of uniform length poly_order + 1
        for idx, func_values in enumerate(values):
            if len(func_values) != poly_order + 1:
                raise ValueError(
                    f"Non-uniform polynomial data at entry {counter}, component {idx}: "
                    f"len(values[{idx}]) = {len(func_values)} but expected {poly_order + 1}"
                )

    check_str(f'n={n}, poly_order={poly_order}')        
    return True


# In[ ]:


def check_str(str_in):
    if VERBOSE:
        print(f"    Structure check passed: {str_in}.")

