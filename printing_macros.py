#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""
This module is used to format prints. 
"""
from sage.rings.real_arb import RealBall
from sage.rings.complex_arb import ComplexBall
from parameters import PRINT_DIGITS, VERBOSE, VERBOSE_COUNTER, V3R8053_C0UN73R


# In[10]:


def print_RBF(ball): 
    """
    PRINTED DATA BY THIS FUNCTION IS NOT RIGOROUS

    Some roundings are done when passing from binary to decimal and rigor is lost.
    In order to avoid this, we always compare the quantities using RBF types.
    """
    # Returns a formatted string "mid ± rad" for RealBall/ComplexBall values.
    if isinstance(ball, (RealBall, ComplexBall)):
        return f"{ball.mid():.{PRINT_DIGITS}g} ± {ball.rad():.{7}g}"
        
    raise TypeError(f"Expected a RealBall/ComplexBall, got {type(ball)}")


# In[ ]:


def print_lemma(verified, lemma_label, bounds, elapsed):
    print("")
    status = "OK" if verified else "FAIL"
    print(f"{lemma_label}: {status}.")
    print(f"Execution time:{print_time(elapsed)}")
    if bounds is not None:
        print("Used bounds were:")
        for (bound, label, side_txt) in bounds:
            if side_txt == 'apprx':
                print(f"    The number {label}: {bound}.")
            else:
                print(f"    The {side_txt} bound {label}: {bound}.")
    print("",flush=True)


# In[ ]:


def print_time(total_seconds):
    three_thou_six_hund = 3600
    sixty = 60
    
    hours = total_seconds // three_thou_six_hund
    minutes_left = total_seconds - hours * three_thou_six_hund
    
    minutes = minutes_left // sixty
    seconds_left = minutes_left - minutes * sixty

    out = ""
    if hours > 0:
        out += f" {hours} hours,"
    if minutes > 0:
        out += f" {minutes} minutes,"
    out += f" {seconds_left} seconds."
    return out


# In[ ]:


def print_iter(counter, text=None, verbose=VERBOSE, verb_count=VERBOSE_COUNTER):
    if not verbose:
        return

    if counter % verb_count != 0:
        return

    indent = "        "
    msg = f"{indent}iter: {counter}"
    if text is not None:
        msg += f", {text}"

    print(msg, flush=True)


# In[ ]:


def print_adapt_inter(tot_iters, tot_ivals, verified_domain, result=None, verbose=VERBOSE, verb_count=V3R8053_C0UN73R):
    if not verbose:
        return

    if tot_iters % verb_count != 0:
        return

    indent = "        "
    msg = f"{indent}Progress: {verified_domain.mid():.4f}, total intervals: {tot_ivals}, total iterations: {tot_iters}"
    if result is not None:
        msg += f", result: {result}"

    print(msg, flush=True)

