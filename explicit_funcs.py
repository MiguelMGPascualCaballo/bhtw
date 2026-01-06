#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from sage.all import *
from printing_macros import *
from classes import FourierRealSeries, Functions_1D
from parameters import *

from collections import deque


# In[ ]:


def equ_tw_eval_symbolic(speed, coeffs, verbose=VERBOSE):
    """
    Computes the residual `xi` of the traveling wave equation (Burgers–Hilbert),
    using symbolic operations on FourierRealSeries.

    In the text:
        xi = c v' + H v + v v'.

    Here:
        - speed is the wave speed c
        - coeffs are cosine Fourier coefficients of v

    Variables:
        Input:
            speed:  (RBF) wave speed c.
            coeffs: (FreeModuleElement / list[RBF]) cosine coefficients of v.
            verbose: (int/bool) verbosity level.
        Output:
            xi: (FourierRealSeries) Fourier series representing xi (typically odd / sine-only).
    """

    #######################################
    # Creates the Fourier series representation v(x) = Σ coeffs[k] cos((k+1)x)
    w = FourierRealSeries.cosine(coeffs)

    #######################################
    # Computes components: v', H v, and v*v'
    dw = w.dx()
    Hw = w.hx()
    quad_term = w * dw
        
    #######################################
    # Computes xi = c*v' + H v + v*v'
    xi = speed * dw + Hw + quad_term

    #######################################
    # Verbose printing
    if verbose:
        printer(f"w is:\n{w}.")
        if verbose >= 2:
            printer(f"dw is:\n{dw}")
            printer(f"Hw is:\n{Hw}")
            printer(f"w * dw  is:\n{quad_term}")
        printer(f"{print_letter('xi')} is:\n{xi}")

    return xi 


# In[ ]:


def w_function_constructor(coeffs, total_derivatives):
    """
    Constructs a Functions_1D object for w from cosine coefficients.

    Variables:
        Input:
            coeffs: (list[RBF]) cosine coefficients of w.
            total_derivatives: (int) number of derivatives to store.
        Output:
            result: (Functions_1D) callable w with derivatives.
    """
    
    #######################################
    # We define 'w' as FourierRealSeries
    w = FourierRealSeries.cosine(coeffs)

    #######################################
    # Build a Functions_1D wrapper and precompute derivatives.
    result = Functions_1D.from_FourierRealSeries(w, total_derivatives)
    result.name = 'w'
    return result


# In[ ]:


def beta_function_constructor(c, w_func):
    """
    Constructs beta(x) = 1 / (c + w(x)).

    Variables:
        Input:
            c:      (RBF) speed parameter.
            w_func: (Functions_1D) w(x) with derivatives.
        Output:
            result: (Functions_1D) beta(x) with derivatives.
    """
    
    denominator =  c + w_func
    result = denominator.inverse(print_letter('beta'))
    return result


# In[ ]:


def betamod_constructor(beta_func, order=4):
    """
    Constructs:
        betamod(x) = (pi - x) * beta(x)^2.

    Implementation detail:
        - aux_betamod(x) = pi - x, with derivatives up to `order`.
        - then betamod = aux_betamod * beta * beta.

    Variables:
        Input:
            beta_func: (Functions_1D) beta(x) with derivatives.
            order: (int) highest derivative order to define for aux_betamod.
        Output:
            result: (Functions_1D) betamod(x) with derivatives.
    """

    #######################################
    # Define aux(x) = pi - x
    aux_betamod = Functions_1D.constant(PI, order) - Functions_1D.identity(tot_ders=order)
    
    #######################################
    # betamod = (pi-x) * beta^2
    result = aux_betamod * beta_func**2
 
    return result


# In[ ]:


def betamod2_constructor(speed, coeffs):
    """
    Intended to build an expression of the form:

        beta^4 + (beta')^2   =  1/(c+w)^4 + (w')^2/(c+w)^4

    Notes:
        - This function currently mixes "speed" into the Fourier series in a way
          that likely does NOT match the intended meaning of w.
        - There are also debug prints left in the function.

    Variables:
        Input:
            speed:  (RBF) wave speed c.
            coeffs: (list[RBF]) cosine coefficients for w.
        Output:
            result: (Functions_1D) enclosure function for the expression above.
    """

    #######################################
    # WARNING:
    # FourierRealSeries(speed, coeffs) means the "mean" term is speed.
    # That would represent (c + w(x)), not w(x).
    ww = FourierRealSeries(speed, coeffs)
    dw = ww.dx()

    #######################################
    # Build (c+w)^4 and (dw)^2
    ww_sq = ww * ww
    dw_sq = dw * dw
    ww_fo = ww_sq * ww_sq
    

    #######################################
    # Numerator: 1 + (dw)^2
    aux_nume = dw_sq + FourierRealSeries.constant(RBF(1))
    
    nume = Functions_1D.from_FourierRealSeries(aux_nume, 4)
    deno = Functions_1D.from_FourierRealSeries(ww_fo, 4)

    #######################################
    # result = (1 + (dw)^2) * (c+w)^(-4)
    result = nume * deno.inverse()

    return result


# In[ ]:


def iota_L1_compute_adaptive(
    integrand,   # this is diota = beta'(x)^2 + beta(x)^4
    beta_sq,     # beta(x)^2
    abs_tol,
    max_iterations=MAX_ITERATIONS,
    verbose=VERBOSE,
    print_each=2**8
):
    """
    Adaptive enclosure for an integral of the form:
        ∫_0^pi beta(x)^2 * gap(x) dx,
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
            beta_sq:   (Functions_1D) beta^2 with derivatives up to order 4.
            abs_tol: (RBF) tolerances scaled by interval radius.
            max_iterations: (int) loop safety bound.
            verbose: (int/bool) verbosity.
        Output:
            result: (RBF) enclosure of the integral.
    """

    #######################################
    # Stack of intervals (LIFO via .pop())
    # We start with [0,pi], and the code processes rightmost parts first
    domain = [ZERO, PI]
    intervals = deque([domain])

    # gap at the RIGHT endpoint of the next interval to process
    # (initially gap(pi)=0)
    last_value = ZERO            # gap at the RIGHT endpoint of the next interval to process

    result = ZERO
    
    current_iterations = 0
    len_domain = RBF(domain[1] - domain[0])
    verified_domain = ZERO

    #######################################
    # order-2 GL remainder uses 4th derivative, so we need n=4 for Leibniz d^4 product
    n = 4
    iota_fat = vector(RBF, n + 1)

    # shorthand for integrand and its 4th derivative enclosure
    g0 = integrand.derivatives[0]
    g4 = integrand.derivatives[4]

    iters = 0
    while intervals:
        iters += 1

        # Optional progress output
        if verbose and (current_iterations % print_each == 0):
            print(f"Total iterations: {current_iterations}, result = {print_RBF(result)}, progress: {verified_domain.mid():.4f}.")
            
        if iters > max_iterations:
            raise RuntimeError("MAX_ITERATIONS reached (no fallback).")
        
        x = intervals.pop()
        current_iterations += 1

        #######################################
        # Interval geometry
        [aa, bb] = x
        
        x0 = (bb + aa) * ONE_DIV_2
        rr = (bb - aa) * ONE_DIV_2

        # 2-pt GL nodes in [a,b]
        x1 = x0 - rr * SQRT_1_DIV_3
        x2 = x0 + rr * SQRT_1_DIV_3

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
        r_ax1 = (x1-aa) * ONE_DIV_2
        r_x12 = (x2-x1) * ONE_DIV_2
        r_x2b = (bb-x2) * ONE_DIV_2
        
        c_ax1 = (x1+aa) * ONE_DIV_2
        c_x12 = (x2+x1) * ONE_DIV_2
        c_x2b = (bb+x2) * ONE_DIV_2

        #######################################
        # GL2 + remainder on each subinterval to enclose the integrals
        int_ax1 = (
            r_ax1 * (g0(c_ax1 + r_ax1*SQRT_1_DIV_3) + g0(c_ax1 - r_ax1*SQRT_1_DIV_3))
            + (r_ax1**5) * g4(RBF([aa, x1])) * ONE_DIV_135
        )
        int_x12 = (
            r_x12 * (g0(c_x12 + r_x12*SQRT_1_DIV_3) + g0(c_x12 - r_x12*SQRT_1_DIV_3))
            + (r_x12**5) * g4(RBF([x1, x2])) * ONE_DIV_135
        )
        int_x2b = (
            r_x2b * (g0(c_x2b + r_x2b*SQRT_1_DIV_3) + g0(c_x2b - r_x2b*SQRT_1_DIV_3))
            + (r_x2b**5) * g4(RBF([x2, bb])) * ONE_DIV_135
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
        # Local GL2 enclosure for ∫_a^b beta_sq(t)*gap(t) dt
        # using nodes x1,x2 and 4th derivative via Leibniz.
        ############################################################

        betaiota_x1 = beta_sq(x1) * iota_th1
        betaiota_x2 = beta_sq(x2) * iota_th2

        beta_fat = [beta_sq.derivatives[idx](RBF(x)) for idx in range(n + 1)]

        # d^4(beta_sq * gap) = Σ_{k=0}^4 C(4,k) beta^(k) * gap^(4-k)
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


# In[ ]:


def palipoly_coeffs_constructor(speed, coeffs):
    """
    Constructs coefficients of a palindromic polynomial encoding:
        P(z) = speed + Σ_{j=1}^rk coeffs[j-1] * (z^j + z^{-j})/2

    This is a common conversion from cosine series to a Laurent polynomial
    in z = e^{ix}.

    Variables:
        Input:
            speed:  (RBF) constant term.
            coeffs: (list[RBF]) cosine coefficients.
        Output:
            polynom: (vector[RBF]) palindromic coefficient vector of length 2*rk+1.
    """
    rk = len(coeffs)
    polynom = zero_vector(RBF,2*rk+1)
    polynom[rk] = speed;
    
    for jj,coef in enumerate(coeffs,start=1):
        polynom[rk-jj] = coef * ONE_DIV_2
        polynom[rk+jj] = coef * ONE_DIV_2

    return polynom


# In[ ]:


def polynom_eval(x, coefficients,idx_der=0):
    """
    Evaluates a polynomial (or its derivative) at x.

    Given coefficients representing:
        P(x) = Σ_{k=0}^{N} coefficients[k] * x^k,
    returns:
        P^{(idx_der)}(x).

    Variables:
        Input:
            x:            (CBF/RBF) evaluation point.
            coefficients: (list/vector) polynomial coefficients.
            idx_der:      (int) derivative order.
        Output:
            result: (CBF) value of the derivative at x.
    """
    
    result = CBF(0)
    for idx_coef, coef in enumerate(coefficients):
        
        if idx_coef>=idx_der:
            #######################################
            # Compute k*(k-1)*...*(k-idx_der+1)
            aux_term = RBF(1)
            
            for jj in range(idx_der):
                aux_term *= idx_coef-jj

            term = aux_term * coef * x**(idx_coef-idx_der)
            result += term

    return result


# In[ ]:


def root_polynomial_is_enclosed(root,polynom,rad,big_bounds):
    """
    Sufficient conditions to verify that a polynomial root is enclosed in a ball.

    Uses Taylor-type bounds with:
        M0 = |P(root)|
        M1 = lower bound on |P'(root)|

    and additional remainder bounds `bound1`, `bound2` built from higher derivative bounds.

    Variables:
        Input:
            root:      (CBF/RBF) approximate root center.
            polynom:   (list/vector) polynomial coefficients.
            rad:       (RBF) candidate enclosure radius.
            big_bounds:(iterable) bounds for higher derivatives (index starts at 1 here).
        Output:
            verified: (bool) True if the enclosure conditions hold.
    """

    #######################################
    # Remainder bounds built from derivatives >= 3
    bound1 = ZERO
    bound2 = ZERO

    for kk,val in enumerate(big_bounds,start=1):
        if kk >= 3:
            bound1 += rad**kk / RBF(factorial(kk)) * val.above_abs()
            bound2 += rad**(kk-1) / RBF(factorial(kk-1)) * val.above_abs()

    #######################################
    # M0 = |P(root)| upper bound, M1 = |P'(root)| lower bound
    M0 = polynom_eval(root,polynom).above_abs()
    M1 = polynom_eval(root,polynom,idx_der=1).below_abs()

    #######################################
    # Enclosure conditions
    cond1 = (M0 + bound1) < rad * M1
    cond2 = bound2 < M1

    return cond1 and cond2
    

