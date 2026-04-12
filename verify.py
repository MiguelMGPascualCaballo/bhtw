#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from sage.rings.real_arb import RealBall
from parameters import RBF, CBF, VERBOSE
from printing_macros import print_RBF, print_check


# In[ ]:





# In[ ]:


def o(out, verbose=VERBOSE):
    """Verifies that out is RBF. In addition, prints out (non-rigorous) if VERBOSE"""
    
    #######################################
    # Verifies that out is a RealBall
    if not isinstance(out, RealBall):
        raise TypeError(f"out must be of type RealBall (RBF) but it is {type(out)}")

    #######################################
    # Verifies that no precision was lost during the computation
    if out.parent() != RBF:
        raise TypeError(f"Precision was lost during the process, out is: {out.parent()}")

    #######################################
    # prints if VERBOSE and returns True (non-rigorous)
    return print_check(f'{print_RBF(out)}', verbose=verbose)


# In[ ]:


def oo(v, n):
    """Verifies that v is an iterable of RBF with n entries"""

    ###################################
    # check that v has n entries
    if len(v) != n:
        raise ValueError(
            f"Length mismatch: v has {len(v)} entries, should be {n}"
        )

    ###################################
    # check that the entries of v are RBF (no printing here)
    for out in v:
        assert o(out, verbose=0)

    ###################################
    # prints if VERBOSE and returns True    
    return print_check(f'iterable of {n} RBF elements')


# In[ ]:


def square_mat_order(M, ord_mat, name_mat='M'):
    """Verifies that M is a square matrix of order ord_mat"""
    
    if M.nrows() != ord_mat or M.ncols() != ord_mat:
        raise ValueError(
            f"{name_mat} must be {ord_mat}x{ord_mat} but is {M.nrows()}x{M.ncols()}."
        )

    ###################################
    # prints if VERBOSE and returns True    
    return print_check(f'{name_mat} is a square matrix of order {ord_mat}')


# In[ ]:


def ball_vector(v):
    """Verifies that v is a Sage vector of RBF or CBF"""
    
    if not hasattr(v, "base_ring"):
        raise TypeError("Input must be a Sage vector")

    base = v.base_ring()
    if base is not RBF and base is not CBF:
        raise TypeError(
            f"Vector base ring must be RBF or CBF, got {base}"
        )
        
    ###################################
    # prints if VERBOSE and returns True
    return print_check(f'vector with base {base}')


# In[ ]:


def data_block(values, n):
    """Verifies that the structure of `ODE approximation data` is correct"""
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

        J, vals = blk

        ###################################
        # vals must have length n
        if len(vals) != n:
            raise ValueError(
                f"Malformed vals entry at index {counter}: "
                f"len(vals) = {len(vals)} but expected {n}"
            )
    
        ###################################
        # Each vals[idx] must be a list of uniform length poly_order + 1
        for idx, func_values in enumerate(vals):
            if len(func_values) != poly_order + 1:
                raise ValueError(
                    f"Non-uniform polynomial data at entry {counter}, component {idx}: "
                    f"len(vals[{idx}]) = {len(func_values)} but expected {poly_order + 1}"
                )

    ###################################
    # prints if VERBOSE and returns True
    return print_check(f'n={n}, poly_order={poly_order}')


# In[ ]:




