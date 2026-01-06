#!/usr/bin/env python
# coding: utf-8

# In[25]:


from sage.all import *
from sage.rings.real_arb import RealBall # Needed to use isinstance

from collections import deque

from methods import *
from parameters import *
from printing_macros import *


# In[26]:


class verify:
    """
    This class provides lightweight sanity checks used throughout the proofs.

    The goal is twofold:
        1) Verify that lemma inputs are given in the expected rigorous numeric type (RBF / RealBall).
        2) Verify that lemma outputs are also RBF and that no precision/type downgrade happened
           during the computation.

    The class is intended to be used as:
        verify.INPUT(speed, coeffs)
        ...
        verify.OUTPUT(result)
    """
    @staticmethod
    def INPUT(c_in, v_in):
        """
        Verifies that the input data of a lemma (for a traveling wave) is correct.

        Specifically:
            - `c_in` must be a RealBall (RBF).
            - `v_in` must be a Python list.
            - every element of `v_in` must be a RealBall (RBF).

        Variables:
            Input:
                c_in: (RBF / RealBall) speed (or constant) of the traveling wave.
                v_in: (list) list of cosine Fourier coefficients, each an RBF.
            Output:
                verified: (bool) True if the verification passes (otherwise raises TypeError).
        """

        #######################################
        # Verifies that c_in is a RealBall (RBF)
        if not isinstance(c_in, RealBall):
            raise TypeError(f"{c_in} must be of type RealBall (RBF)")
    
        #######################################
        # Verifies that v_in is a list
        if not isinstance(v_in, list):
            raise TypeError(f"{v_in} must be a list")
    
        #######################################
        # Verifies that all elements of v_in are RealBall (RBF)
        for ii, element in enumerate(v_in):
            if not isinstance(element, RealBall):
                raise TypeError(f"Element {ii} in v_in is not a RealBall (RBF).")

        #######################################
        # Logging / progress message
        printer_initial(f"Input verified, computing for speed: c={print_RBF(c_in)}.")
        
        return True

    @staticmethod
    def OUTPUT(out):
        """
        Verifies that the output of a lemma is correct.

        Specifically:
            - `out` must be a RealBall (RBF).
            - `out.parent()` must be exactly `RBF`, to ensure no precision/type downgrade happened.

        Variables:
            Input:
                out: (RBF / RealBall) result produced by the lemma.
            Output:
                verified: (bool) True if the verification passes (otherwise raises TypeError).
        """
        
        #######################################
        # Verifies that out is a RealBall (RBF)
        if not isinstance(out, RealBall):
            raise TypeError(f"OUTPUT must be of type RealBall (RBF) but it is {type(out)}")
    
        #######################################
        # Verifies that no precision was lost during the computation
        # (i.e., the parent ring/type is still RBF)
        if out.parent() != RBF:
            raise TypeError(f"Precision was lost during the process, out is: {out.parent()}")

        #######################################
        # Logging / progress message
        printer_initial(f"Output verified: {print_RBF(out)}")
        
        return True


# In[ ]:


def sup_bound_prev_esti_adapt_1D(
    domain, func, guess_bound, method,
    max_iterations=MAX_ITERATIONS,
    max_subintervals=MAX_SUBINTERVALS,
    verbose=VERBOSE,
    print_each=PRINT_EACH
):
    """
    Verifies an a priori supremum bound using adaptive domain subdivision.

    The function checks whether the strict inequality
        sup_{x in domain} |func(x)| < guess_bound
    holds, using interval/ball arithmetic. A user-provided local enclosure method
    `method` must satisfy: for every subinterval I ⊂ domain,
        y = method(I, func)
    returns an enclosure y such that func(I) ⊂ y.

    Adaptive strategy:
      - Start with the full domain.
      - For a subinterval I, compute y = method(I, func).
        * If y.abs() < guess_bound, the bound is verified on I.
        * If y.abs() > guess_bound, the bound is violated on I and the function returns False.
        * Otherwise, the test is inconclusive and I is split into two halves.

    Stopping outcomes:
        returns True
            if the strict bound is verified on the entire domain.
        returns False
            if a subinterval is proven to violate the bound.
        raises RuntimeError("MAX_ITERATIONS_REACHED")
            if the total number of calls to `method` exceeds `max_iterations`.
        raises RuntimeError("MAX_INTERVALS_REACHED")
            if the number of active subintervals exceeds `max_subintervals`.

    Notes:
      - This routine verifies a strict inequality (<). If the true supremum equals
        `guess_bound` exactly, the algorithm may fail to conclude and may hit a limit.
      - Progress messages are printed only if `verbose` is enabled.

    Parameters:
        domain (RBF interval):
            One-dimensional domain on which the supremum is tested (e.g. RBF([0, PI])).
        func (callable):
            Function to bound; must be compatible with `method`.
        guess_bound (RBF):
            Proposed upper bound for sup_{x in domain} |func(x)|.
        method (callable):
            Local enclosure routine y = method(I, func) returning an enclosure of func(I).

        max_iterations (int, optional):
            Maximum number of method evaluations allowed before aborting.
        max_subintervals (int, optional):
            Maximum number of subintervals allowed to be stored during subdivision.
        verbose (bool or int, optional):
            If truthy, prints progress information during the computation.
        print_each (int, optional):
            When `verbose` is enabled, print progress every `print_each` steps.

    Returns:
        bool:
            True iff the algorithm verifies sup_{x in domain} |func(x)| < guess_bound.
    """

    # Stack of subintervals to process (we adaptively refine until every interval is verified)
    intervals = deque([domain])

    # Counters used for progress printing
    current_iterations = 0
    len_domain = RBF(domain[1] - domain[0])
    verified_domain = ZERO
    #################################################
    while intervals:

        # Optional progress output
        if verbose and (len(intervals) % print_each == 0):
            print(f"Total intervals: {len(intervals)}, progress: {verified_domain.mid():.4f}.")
        if verbose and (current_iterations % print_each == 0):
            print(f"Total iterations: {current_iterations}, progress: {verified_domain.mid():.4f}.")

        # Pop a subinterval and compute an enclosure of func on it
        x = intervals.pop()
        y = method(x, func)
        current_iterations += 1

        # Case 1: provably inside the bound on this subinterval
        if y.abs() < guess_bound:
            verified_domain += RBF(x[1] - x[0]) / len_domain

        # Case 2: provably violates the bound on this subinterval
        elif y.abs() > guess_bound:
            return False

        # Case 3: inconclusive -> subdivide and retry on smaller intervals
        else: 
            aa, bb = x
            x0 = (bb + aa) * ONE_DIV_2
            
            xl = [aa, x0]
            xr = [x0, bb]

            intervals.append(xl)
            intervals.append(xr)

        # Safety checks to prevent runaway refinement
        if len(intervals) > max_subintervals:
            raise RuntimeError("MAX_INTERVALS_REACHED")
        if current_iterations > max_iterations:
            raise RuntimeError("MAX_ITERATIONS_REACHED")

    # Optional progress output
    if verbose:
        print(
            f" Total intervals: {len(intervals)}.\n Total iterations: {current_iterations}.\n Progress: {verified_domain.mid():.4f}.\n"
        )
        
    # If we exit the loop, every subinterval has been verified
    return True


# In[ ]:


def integral_adaptive_computation_1D(
    domain, func, method,
    abs_tol=ABS_TOL_1D,
    rel_tol=REL_TOL_1D,
    max_iterations=MAX_ITERATIONS,
    max_subintervals=MAX_SUBINTERVALS,
    verbose=VERBOSE,
    print_each=PRINT_EACH
):
    """
    Computes an enclosure of a one-dimensional quantity using adaptive subdivision.

    This routine adaptively subdivides `domain` and applies a user-provided local
    enclosure method `method` on each subinterval. The returned value is the sum of
    the accepted local enclosures:
        result = Σ y_I,
    where y_I = method(I, func) is the enclosure produced on subinterval I.

    In typical usage, `method` is an integral enclosure routine, so that y_I encloses
    ∫_I func(x) dx and the final result encloses ∫_domain func(x) dx. (More generally,
    the method may compute any additive quantity defined in your `methods` module.)

    Acceptance and refinement:
      - For a subinterval I, compute y = method(I, func).
      - Let r_y = rad(y), r_x = rad(I), and c_y be the center of y.
      - Accept y (add it to the result) if either:
            r_y < abs_tol * r_x                      (absolute criterion), or
            r_y < rel_tol * r_x * |c_y|              (relative criterion).
      - Otherwise split I into two halves and continue.

    Fallback mode:
      - If the number of queued subintervals exceeds `max_subintervals`, or the number
        of method evaluations exceeds `max_iterations`, the algorithm enters fallback mode.
      - In fallback mode, no further subdivision is performed; remaining subintervals are
        evaluated once and accumulated. The result remains a valid enclosure, possibly less sharp.

    Parameters:
        domain (RBF interval):
            One-dimensional domain on which the computation is performed.
        func (callable):
            Function to evaluate; must be compatible with `method`.
        method (callable):
            Local enclosure routine y = method(I, func) returning an enclosure of the
            subinterval contribution.

        abs_tol (RBF or float, optional):
            Absolute tolerance controlling when a subinterval contribution is accepted.
        rel_tol (RBF or float, optional):
            Relative tolerance controlling when a subinterval contribution is accepted.
        max_iterations (int, optional):
            Maximum number of calls to `method` allowed before entering fallback mode.
        max_subintervals (int, optional):
            Maximum number of queued subintervals allowed before entering fallback mode.
        verbose (bool or int, optional):
            If truthy, prints progress information during the computation.
        print_each (int, optional):
            When `verbose` is enabled, print progress every `print_each` steps.

    Returns:
        RBF:
            Enclosure of the computed additive quantity over the full domain
            (typically an enclosure of ∫_domain func(x) dx).
    """
    
    # Queue of subintervals to process (we adaptively refine until every interval is verified)
    intervals = deque([domain])

    # Running enclosure of the integral (additive aggregation)
    result = ZERO

    # Counter for method evaluations (used for progress printing and safety limits)
    current_iterations = 0
    len_domain = RBF(domain[1] - domain[0])
    verified_domain = ZERO

    # Once this is True, we stop subdividing and simply accumulate remaining intervals.
    fallback_mode = False
    
    #################################################
    while intervals:

        # Optional progress output
        if verbose and (len(intervals) % print_each == 0):
            print(f"Total intervals: {len(intervals)}, result = {print_RBF(result)}, progress: {verified_domain.mid():.4f}.")
        if verbose and (current_iterations % print_each == 0):
            print(f"Total iterations: {current_iterations}, result = {print_RBF(result)}, progress: {verified_domain.mid():.4f}.")

        # Take the next interval
        x = intervals.pop()
        y = method(x, func)
        current_iterations += 1
        
        # Radii used in acceptance tests
        r_y = RBF(y.rad())          # enclosure radius of local contribution
        c_y_abs = abs(y.squash())   # |center(y)| for relative scaling
        d_x = x[1] - x[0]           # length of the interval

        # In fallback mode we accept everything without further subdivision
        if fallback_mode:
            result += y
            verified_domain += d_x / len_domain
            continue
        
        # Decide whether this interval is resolved to the requested tolerance
        crit_abs = r_y < abs_tol * d_x
        crit_rel = r_y < rel_tol * d_x * c_y_abs

        # Accept if the enclosure is sufficiently tight
        if crit_abs or crit_rel:
            result += y
            verified_domain += d_x / len_domain
            continue
            
        # Otherwise split and enqueue children (to the left)
        aa, bb = x
        x0 = (bb + aa) * ONE_DIV_2
        xl = [aa, x0]
        xr = [x0, bb]

        intervals.appendleft(xl)
        intervals.appendleft(xr)

        # Switch to fallback mode (no more splitting) if too many intervals or iterations
        if len(intervals) > max_subintervals:
            fallback_mode = True
            if verbose:
                print("MAX_SUBINTERVALS reached — entering fallback mode.")
        elif current_iterations > max_iterations:
            fallback_mode = True
            if verbose:
                print("MAX_ITERATIONS reached — entering fallback mode.")

    # Optional progress output
    if verbose:
        print(
            f" Total intervals: {len(intervals)}.\n Total iterations: {current_iterations}.\n Progress: {verified_domain.mid():.4f}.\n"
        )
    
    return result

