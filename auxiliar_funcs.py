#!/usr/bin/env python
# coding: utf-8

# In[25]:


from sage.all import *
from sage.rings.real_arb import RealBall # Needed to use isinstance
from sage.modules.free_module_element import vector

from collections import deque

from parameters import (
    RBF, CBF,
    ZERO, TWO, PI, ONE_DIV_2, SQRT_1_DIV_3, ONE_DIV_135,
    VERBOSE, V3R8053_C0UN73R,
    ABS_TOL, REL_TOL, MAX_ITERATIONS, MAX_SUBINTERVALS,
)
from printing_macros import print_RBF, print_adapt_inter
from load_data import interval_from_J
import verify


# In[ ]:


def norm_vector_sq(v):
    """Squared norm of a Sage RBF or CBF vector: sum |v_k|^2, returned as RBF."""
    return prod_vector(v, v).real()


# In[ ]:


def prod_vector(u, v):
    """Hermitian inner product <u, v> = sum u_i * conj(v_i), returned as RBF (real) or CBF."""
    verify.ball_vector(u)
    verify.ball_vector(v)
    
    assert len(u) == len(v)

    base = v.base_ring()

    if base is RBF:
        return sum(a * b for (a, b) in zip(u, v))

    return sum(a * b.conjugate() for (a, b) in zip(u, v))


# In[ ]:


def build_CS_for_jj(jj, Rk, x0, rr):
    """
    Computes the Taylor-like moment lists C and S of orders 0..Rk for the interval
    [x0-rr, x0+rr] and frequency jj, defined by:
        C[l] = ∫_{x0-rr}^{x0+rr} (x-x0)^l cos(jj*x) dx
        S[l] = ∫_{x0-rr}^{x0+rr} (x-x0)^l sin(jj*x) dx

    Computed via recurrence in l. Negative jj is handled by symmetry (S flips sign).
    Returns (C, S), each a list of RBF balls of length Rk+1, or ([], []) if Rk < 0.
    """
    
    if Rk < 0:
        return [], []

    signS = 1
    if jj < 0:
        jj = -jj
        signS = -1

    invj = 1 / RBF(jj)

    aa = x0 - rr
    bb = x0 + rr

    C = [None] * (Rk + 1)
    S = [None] * (Rk + 1)

    if jj == 0:
        for ll in range(Rk + 1):
            if ll % 2 == 0:
                C[ll] = 2 * (rr ** (ll + 1)) / RBF(ll + 1)
            else:
                C[ll] = 0
            S[ll] = 0
        return C, S 

    # k=0
    C[0] = (sin(jj * bb) - sin(jj * aa)) * invj
    S[0] = (cos(jj * aa) - cos(jj * bb)) * invj

    # Recurrence
    for ll in range(1, Rk + 1):
        s_term = rr**ll * invj * (sin(jj * bb) - ((-1) ** ll) * sin(jj * aa))
        c_term = rr**ll * invj * (((-1) ** ll) * cos(jj * aa) - cos(jj * bb))

        C[ll] = s_term - ll * invj * S[ll - 1]
        S[ll] = c_term + ll * invj * C[ll - 1]

    # Change sign of S if jj was negative
    if signS == -1:
        S = [ -v for v in S ]

    return C, S


# In[ ]:


def precompute_expo_vals(J, n_min, total, max_Nk):
    """
    Precomputes the CBF matrix expo_vals of shape (total-n_min) x max_Nk, where
        expo_vals[jj-n_min, kk] = C[kk] - i*S[kk]
    with C, S the moment lists from build_CS_for_jj for the interval J and frequency jj.
    Used to assemble Fourier-type operators efficiently.
    """
    _, _, x0, rr = interval_from_J(J)
    
    expo_vals = matrix(CBF, total - n_min, max_Nk)
    for jj in range(n_min, total):
        idx_jj = jj-n_min 
        C, S = build_CS_for_jj(jj, max_Nk, x0=x0, rr=rr)
        for kk in range(max_Nk):
            expo_vals[idx_jj, kk] = C[kk] - I*S[kk]

    return expo_vals



# In[ ]:


def s_from_24(V):
    """
    Measures how far V is from unitary: returns s = max_{j<=k} |G[j,k]|  where  G = V^H V - I.
    Raises if s exceeds 1/(8*n), which is the threshold required by Lemma 2.4 of:
        `ANY THREE EIGENVALUES DO NOT DETERMINE A TRIANGLE`
    """
    ord_mat = V.nrows()
    assert ord_mat == V.ncols()

    threshold = 1 / RBF(8 * ord_mat)

    # Compute G = V* V - I in one matrix multiplication
    base = V.base_ring()
    if base is CBF:
        Vh = V.conjugate().transpose()
    else:
        Vh = V.transpose()
    G = Vh * V - matrix.identity(CBF, ord_mat)

    # Find the max |entry| over upper triangle (including diagonal)
    s = ZERO
    for jj in range(ord_mat):
        for kk in range(jj + 1):
            a = G[kk, jj].abs()
            if a > s:
                s = a
                if not s < threshold:
                    raise Exception("Assumption not satisfied.")

    return RBF(s.upper())


# In[ ]:


def S_from_ti_S_for_Gersh(ti_S, M, V):
    """
    Given ti_S = V^H (M* M) V and the approximate eigenvector matrix V,
    returns an inflated matrix S such that every eigenvalue of M^* M
    lies inside the union of Gershgorin discs of S.

    The inflation uses Weyl's perturbation bound: each entry M[j,k] is widened by
        s1*(||MV col_j|| + ||MV col_k||) + s2*||M||_F,
    where s = s_from_24(V), s1 = sqrt(3s), s2 = 4s.
    """
    ord_mat = ti_S.ncols()
    assert ord_mat == ti_S.nrows()
    assert ord_mat == M.nrows()
    assert ord_mat == M.ncols()
    assert ord_mat == V.nrows()
    assert ord_mat == V.ncols()

    s = s_from_24(V)
    
    s1 = (3*s)**ONE_DIV_2
    s2 = 4*s
    
    MV = M * V
    MV_norms = [col.norm() for col in MV.columns()] # function from explicit_funcs
    MF = fro_norm(M)

    base_error = s2 * MF  # scalar, same for all entries
    S = matrix(CBF, ord_mat)
    for jj, Mvj in enumerate(MV_norms):
        row_base = s1 * Mvj + base_error
        for kk, Mvk in enumerate(MV_norms):
            term = ti_S[jj, kk]
            rr = row_base + s1 * Mvk
            S[jj, kk] = term.add_error(rr)

    return S


# In[ ]:


def gersh_smallest_eig(M):
    """
    Gershgorin lower bound for the smallest eigenvalue of the Hermitian matrix M:
        min_j ( Re(M[j,j]) - sum_{k!=j} |M[j,k]| ).
    """
    ord_mat = M.ncols()
    assert ord_mat == M.nrows()
    
    vals = vector(RBF, ord_mat)
    for jj in range(ord_mat):
        val = M[jj, jj]
        for kk in range(ord_mat):
            if kk == jj:
                continue
            val -= M[jj, kk].abs()
        vals[jj] = val.real()
        
    return min(vals)


# In[ ]:


def gersh_gap_one_two(M):
    """
    Gershgorin bounds for the spectral gap between the smallest and second-smallest
    eigenvalue of the Hermitian matrix M, assumed sorted in DECREASING diagonal order.   <---- IMPORTANT

    Returns (val1, val2) where:
        val1 = upper bound for the smallest eigenvalue  (last Gershgorin disc upper edge)
        val2 = lower bound for the second eigenvalue    (min of remaining Gershgorin lower edges)

    A spectral gap is rigorously certified when val2 > val1 is proven.
    """
    ord_mat = M.ncols()
    assert ord_mat == M.nrows()
    
    val1 =  M[-1, -1]
    for kk in range(ord_mat-1):
        val1 += M[-1, kk].abs()
    val1 = val1.real()

    vals2 = vector(RBF, ord_mat-1)
    for jj in range(ord_mat-1):
        val = M[jj, jj]
        for kk in range(ord_mat):
            if kk == jj:
                continue
            val -= M[jj, kk].abs()
        vals2[jj] = val.real()

    return val1, min(vals2)


# In[ ]:


def fro_norm(A):
    """Frobenius norm of A (real or complex CBF matrix): ( sum_{i,j} |A_ij|^2 )^{1/2}."""
    
    result = ZERO
    for idx_row in range(A.nrows()):
        for idx_col in range(A.ncols()):
            z = A[idx_row, idx_col] 
            to_add = (z * z.conjugate()).real()
            result += to_add
            
    return result ** ONE_DIV_2


# In[ ]:


def max_prev_esti(
    domain, func, guess_bound, method,
    max_iterations=MAX_ITERATIONS,
    max_subintervals=MAX_SUBINTERVALS,
    verbose=VERBOSE,
    verb_count=V3R8053_C0UN73R
):
    """
    Verifies an a priori supremum bound using adaptive domain subdivision.

    The function checks whether the strict inequality
        sup_{x in domain} func(x) < guess_bound
    holds, using interval/ball arithmetic. A user-provided local enclosure method
    `method` must satisfy: for every subinterval I ⊂ domain,
        y = method(I, func)
    returns an enclosure y such that func(I) ⊂ y.

    Adaptive strategy:
      - Start with the full domain.
      - For a subinterval I, compute y = method(I, func).
        * If y.abs < guess_bound, the bound is verified on I.
        * If y.abs > guess_bound, the bound is violated on I and the function returns False.
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
            Proposed upper bound for sup_{x in domain} func(x).
        method (callable):
            Local enclosure routine y = method(I, func) returning an enclosure of func(I).

        max_iterations (int, optional):
            Maximum number of method evaluations allowed before aborting.
        max_subintervals (int, optional):
            Maximum number of subintervals allowed to be stored during subdivision.
        verbose (bool or int, optional):
            If truthy, prints progress information during the computation.
        verb_count (int, optional):
            When `verbose` is enabled, print progress every `print_each` steps.

    Returns:
        bool:
            True iff the algorithm verifies sup_{x in domain} func(x) < guess_bound.
    """
    half = ONE_DIV_2
    # Stack of subintervals to process (we adaptively refine until every interval is verified)
    intervals = deque([domain])

    # Counters used for progress printing
    
    len_domain = RBF(domain[1] - domain[0])
    verified_domain = ZERO
    current_iterations = 0
    while intervals:
        current_iterations += 1

        # Optional progress output
        print_adapt_inter(current_iterations, len(intervals), verified_domain, verbose=verbose, verb_count=verb_count)

        # Pop a subinterval and compute an enclosure of func on it
        x = intervals.pop()
        y = method(x, func)
        
        # Case 1: provably inside the bound on this subinterval
        if y < guess_bound:
            verified_domain += RBF(x[1] - x[0]) / len_domain

        # Case 2: provably violates the bound on this subinterval
        elif y > guess_bound:
            return False

        # Case 3: inconclusive -> subdivide and retry on smaller intervals
        else: 
            aa, bb = x
            x0 = (bb + aa) * half
            
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


def integ_adaptive_1D(
    domain, func, method,
    abs_tol=ABS_TOL,
    rel_tol=REL_TOL,
    max_iterations=MAX_ITERATIONS,
    max_subintervals=MAX_SUBINTERVALS,
    verbose=VERBOSE,
    verb_count=V3R8053_C0UN73R
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
        verb_count (int, optional):
            When `verbose` is enabled, print progress every `print_each` steps.

    Returns:
        RBF:
            Enclosure of the computed additive quantity over the full domain
            (typically an enclosure of ∫_domain func(x) dx).
    """
    
    half = ONE_DIV_2
    
    # Queue of subintervals to process (we adaptively refine until every interval is verified)
    intervals = deque([domain])
    len_domain = RBF(domain[1] - domain[0])
    
    fallback_mode = False    # Once this is True, we stop subdividing and simply accumulate remaining intervals.
    current_iterations = 0   # Counter for method evaluations (used for progress printing and safety limits)
    verified_domain = ZERO
    result = ZERO            # Running enclosure of the integral (additive aggregation)
    #################################################
    while intervals:

        # Optional progress output
        print_adapt_inter(current_iterations, len(intervals), verified_domain, result=print_RBF(result), verbose=verbose, verb_count=verb_count)

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
        x0 = (bb + aa) * half
        xl = [aa, x0]
        xr = [x0, bb]

        intervals.appendleft(xl)
        intervals.appendleft(xr)

        # Switch to fallback mode (no more splitting) if too many intervals or iterations
        if len(intervals) > max_subintervals:
            fallback_mode = True
            print("MAX_SUBINTERVALS reached — entering fallback mode.")
        elif current_iterations > max_iterations:
            fallback_mode = True
            print("MAX_ITERATIONS reached — entering fallback mode.")

    # Optional progress output
    if verbose:
        print(
            f" Total intervals: {len(intervals)}.\n Total iterations: {current_iterations}.\n Progress: {verified_domain.mid():.4f}.\n"
        )
    
    return result


# In[ ]:


def triangular_integral(
    integrand,   # this is diota = beta'(x)^2 + beta(x)^4
    mult_term,   # beta(x)^2
    abs_tol=ABS_TOL,
    max_iterations=MAX_ITERATIONS,
    verbose=VERBOSE,
    verb_count=V3R8053_C0UN73R
):
    """
    Adaptive enclosure for an integral of the form:
        ∫_0^pi mult_term(x) * gap(x) dx,
    where:
        gap(x) = ∫_x^pi integrand(t) dt
    so:
        gap(pi)=0 and gap'(x) = -integrand(x).

    The algorithm:
        - Maintains a stack of intervals in [0,pi] and processes from right to left.
        - On each interval [a,b], constructs enclosures of gap at quadrature nodes.
        - Uses 2-point Gauss–Legendre + 4th-derivative remainder (order-2 GL with remainder).
        - If local error enclosure is small enough, accept; otherwise bisect.

    Variables:
        Input:
            integrand: (Functions_1D) must provide derivatives up to order 4.
            mult_term:   (Functions_1D) beta^2 with derivatives up to order 4.
            abs_tol: (RBF) tolerances scaled by interval radius.
            max_iterations: (int) loop safety bound.
            verbose (bool or int, optional): If truthy, prints progress information during the computation.
            verb_count (int, optional): When verbose is enabled, print progress every `print_each` steps.
        Output:
            result: (RBF) enclosure of the integral.
    """

    half = ONE_DIV_2
    w3 = SQRT_1_DIV_3
    o135 = ONE_DIV_135
    
    #######################################
    # Stack of intervals (LIFO via .pop())
    # We start with [0,pi], and the code processes rightmost parts first
    domain = [ZERO, PI]
    intervals = deque([domain])

    # gap at the RIGHT endpoint of the next interval to process
    # (initially gap(pi)=0)
    last_value = ZERO            # gap at the RIGHT endpoint of the next interval to process

    result = ZERO
    
    len_domain = RBF(domain[1] - domain[0])
    verified_domain = ZERO

    #######################################
    # order-2 GL remainder uses 4th derivative, so we need n=4 for Leibniz d^4 product
    n = 4
    iota_fat = vector(RBF, n + 1)

    # shorthand for integrand and its 4th derivative enclosure
    g0 = integrand.derivatives[0]
    g4 = integrand.derivatives[4]

    current_iterations = 0
    while intervals:
        current_iterations += 1

        # Optional progress output
        print_adapt_inter(current_iterations, len(intervals), verified_domain, result=print_RBF(result), verbose=verbose, verb_count=verb_count)

        if current_iterations > max_iterations:
            raise RuntimeError("MAX_ITERATIONS reached (no fallback).")
        
        x = intervals.pop()

        #######################################
        # Interval geometry
        [aa, bb] = x
        
        x0 = (bb + aa) * half
        rr = (bb - aa) * half

        # 2-pt GL nodes in [a,b]
        x1 = x0 - rr * w3
        x2 = x0 + rr * w3

        ############################################################
        # Build gap values on this interval FROM THE RIGHT:
        # ib = gap(b) is known as last_value.
        #
        # gap(x2) = gap(b) + ∫_{x2}^{b} integrand
        # gap(x1) = gap(x2) + ∫_{x1}^{x2} integrand
        # gap(a)  = gap(x1) + ∫_{a}^{x1} integrand
        ############################################################

        #######################################
        # Split [a,b] into 3 subintervals: [a,x1], [x1,x2], [x2,b]
        r_ax1 = (x1-aa) * half
        r_x12 = (x2-x1) * half
        r_x2b = (bb-x2) * half
        
        c_ax1 = (x1+aa) * half
        c_x12 = (x2+x1) * half
        c_x2b = (bb+x2) * half

        #######################################
        # GL2 + remainder on each subinterval to enclose the integrals
        int_ax1 = (
            r_ax1 * (g0(c_ax1 + r_ax1*w3) + g0(c_ax1 - r_ax1*w3))
            + (r_ax1**5) * g4(RBF([aa, x1])) * o135
        )
        int_x12 = (
            r_x12 * (g0(c_x12 + r_x12*w3) + g0(c_x12 - r_x12*w3))
            + (r_x12**5) * g4(RBF([x1, x2])) * o135
        )
        int_x2b = (
            r_x2b * (g0(c_x2b + r_x2b*w3) + g0(c_x2b - r_x2b*w3))
            + (r_x2b**5) * g4(RBF([x2, bb])) * o135
        )

        #######################################
        # Gap values at key points
        ib  = last_value          # gap(b)
        ix2 = ib + int_x2b        # gap(x2)
        ix1 = ix2 + int_x12       # gap(x1)
        ia  = ix1 + int_ax1       # gap(a)
        
        # gap enclosure on [a,b] (gap decreasing): [gap(b), gap(a)]
        ix = RBF([ib, ia])

        iota_th1 = ix1
        iota_th2 = ix2

        #######################################
        # iota_fat[0] is gap([a,b]) enclosure
        iota_fat[0] = ix
        
        # since gap'(x) = -integrand(x), we have
        # gap^(k)(x) = - integrand^(k-1)(x) for k>=1
        for idx_der in range(1, n + 1):
            g = integrand.derivatives[idx_der - 1]
            iota_fat[idx_der] = -g(RBF(x))   # interval enclosure on x

        ############################################################
        # Local GL2 enclosure for ∫_a^b mult_term(t)*gap(t) dt
        # using nodes x1,x2 and 4th derivative via Leibniz.
        ############################################################

        betaiota_x1 = mult_term(x1) * iota_th1
        betaiota_x2 = mult_term(x2) * iota_th2

        beta_fat = [mult_term.derivatives[idx](RBF(x)) for idx in range(n + 1)]

        # d^4(mult_term * gap) = Σ_{k=0}^4 C(4,k) beta^(k) * gap^(4-k)
        aux_vec = [1, 4, 6, 4, 1]
        betaiota_d4 = sum(aux_vec[k] * beta_fat[k] * iota_fat[4 - k] for k in range(5))

        y = rr * (betaiota_x1 + betaiota_x2) + betaiota_d4 * (rr**5) / RBF(135)

        ############################################################
        # Acceptance / refinement
        ############################################################

        c_y_abs = y.squash().abs()
        r_y = RBF(y.rad())

        crit_abs = r_y < abs_tol * rr

        if crit_abs:
            # IMPORTANT: update last_value to gap(a),
            # because next interval (to the left) has right endpoint = a
            last_value = ia
            result += y
            verified_domain += 2*rr / len_domain
            continue

        #######################################
        # Refine: split [a,b] into [a,x0] and [x0,b]
        xl = [aa, x0]
        xr = [x0, bb]

        # push left THEN right so right is processed next (.pop())
        intervals.append(xl)
        intervals.append(xr)

    # Optional progress output
    if verbose:
        print(
            f" Total intervals: {len(intervals)}.\n Total iterations: {current_iterations}.\n Progress: {verified_domain.mid():.4f}.\n"
        )

    return result

