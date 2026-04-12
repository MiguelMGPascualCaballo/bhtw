#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from sage.all import *
from classes import FourierRealSeries, Functions_1D
from parameters import (
    RBF, CBF, 
    ZERO, ONE, TWO, 
    PI, TWOPI, ONE_DIV_2, 
)
from load_data import interval_from_J
import auxiliar_funcs
import verify
from printing_macros import print_iter_J

from collections import deque


# In[ ]:


def equ_tw_eval_symbolic(speed, coeffs):
    """
    Computes the traveling wave residual xi = c*v' + H(v) + v*v'
    as a FourierRealSeries, where v = sum coeffs[k]*cos((k+1)x).
    The result is typically sine-only (odd function).
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

    return xi 


# In[ ]:


def fvap_der_and_Hil(fvap, N):
    """
    Returns the Fourier coeffients of the derivatives and Hilbert transform of fvap.
    Here we have that:
        fvap = [0, 1, ..., N-1, -1, -2, ..., -N]
    """
    def _stab_mult(func, N): 
        aux_pos = [func(   jj) for jj in range(N)]
        aux_neg = [func(-1-jj) for jj in range(N)]
        return vector(CBF, aux_pos + aux_neg)

    def _coeffs_diff(jj):
        return I*RBF(jj)

    def _coeffs_Hilb(jj):
        result = 0
        if jj >= 1:
            result = 1
        elif jj <= -1:
            result = -1
        return -I * result

    def _zip_mult(f, g): 
        assert len(f) == len(g)
        return vector(CBF, [a * b for a, b in zip(f, g)])
        
    der_coeffs = _stab_mult(_coeffs_diff, N)
    Hil_coeffs = _stab_mult(_coeffs_Hilb, N)
    
    fvap_der = _zip_mult(fvap, der_coeffs)
    fvap_Hil = _zip_mult(fvap, Hil_coeffs)

    return fvap_der, fvap_Hil


# In[ ]:


def beta_func_constructor(speed, coeffs, order=4):
    """Constructs beta(x) = 1 / (speed + w(x)) as a Functions_1D, where w = sum coeffs[k]*cos((k+1)x)."""
    ww = FourierRealSeries(speed, coeffs)
    deno = Functions_1D.from_FourierRealSeries(ww, order)
    return deno.inverse()
    


# In[ ]:


def aux_kappa1_func_constructor(speed, coeffs, order=4):
    """Constructs (pi - x) * beta(x)^2 = (pi - x) / (speed + w(x))^2 as a Functions_1D."""

    #######################################
    # Define aux(x) = pi - x
    aux_func = Functions_1D.constant(PI, order) - Functions_1D.identity(tot_ders=order)

    #######################################
    # Define beta(x)^2 = 1 / (c+v(x))^2
    be_sq = beta_sq_func_constructor(speed, coeffs, order)
    
    #######################################
    # betamod = (pi-x) * beta^2
    result = aux_func * be_sq
    return result


# In[ ]:


def betamod_func_sq_constructor(speed, coeffs, order):
    """Constructs (1 + w'(x)^2) / (speed + w(x))^4 as a Functions_1D."""

    #######################################
    ww = FourierRealSeries(speed, coeffs)
    dw = ww.dx()

    #######################################
    # Build (c+w)^4 and (dw)^2
    ww_sq = ww * ww
    dw_sq = dw * dw
    ww_fo = ww_sq * ww_sq

    #######################################
    # Numerator: 1 + (dw)^2
    aux_nume = dw_sq + FourierRealSeries.constant(ONE)
    
    nume = Functions_1D.from_FourierRealSeries(aux_nume, order)
    deno = Functions_1D.from_FourierRealSeries(ww_fo, order)

    #######################################
    # result = (1 + (dw)^2) * (c+w)^(-4)
    return nume * deno.inverse()


# In[ ]:


def beta_dx_func_sq_constructor(speed, coeffs, order):
    """Constructs w'(x)^2 / (speed + w(x))^4 as a Functions_1D."""

    #######################################
    ww = FourierRealSeries(speed, coeffs)
    dw = ww.dx()

    #######################################
    # Build (c+w)^4 and (dw)^2
    ww_sq = ww * ww
    dw_sq = dw * dw
    ww_fo = ww_sq * ww_sq
    
    nume = Functions_1D.from_FourierRealSeries(dw_sq, order)
    deno = Functions_1D.from_FourierRealSeries(ww_fo, order)

    #######################################
    # result = (1 + (dw)^2) * (c+w)^(-4)
    return nume * deno.inverse()


# In[ ]:


def real_gk_funcs_constructor(coeffs, order):
    N = len(coeffs)
    
    func_list = []
    for kk in range(N):
        sine_coeffs = [vkm for vkm in coeffs[kk:]]
        sine_series = FourierRealSeries.sine(sine_coeffs) * RBF(kk)
        func = Functions_1D.from_FourierRealSeries(sine_series, order)
        func_list.append(func)

    return func_list


# In[ ]:


def beta_sq_func_constructor(speed, coeffs, order):
    """Constructs beta(x)^2 = 1 / (speed + w(x))^2 as a Functions_1D."""
    ww = FourierRealSeries(speed, coeffs)
    ww_sq = ww * ww
    return  Functions_1D.from_FourierRealSeries(ww_sq, order).inverse()


# In[ ]:


def betamod_laap_sq_func_constructor(speed, coeffs, laap, the, order):
    """
    Constructs the pair (f_+, f_-) where
        f_±(x) = |(i ± la)*beta(x) - i*the|^2,
    with beta(x) = 1/(speed + w(x)). Returns two Functions_1D objects.
    """
    ww = FourierRealSeries(speed, coeffs)
    ww_sq = ww * ww

    lare = laap.real()
    laim = laap.imag()
    la_norm_sq = lare**2 + laim**2
    
    A_p = la_norm_sq + 1 + 2*laim
    A_m = la_norm_sq + 1 - 2*laim
    B_p = -2*the*(laim + 1)
    B_m = -2*the*(laim - 1)
    C = the**2

    aux_nume_p = A_p + B_p*ww + C*ww_sq
    aux_nume_m = A_m + B_m*ww + C*ww_sq
    
    nume_p = Functions_1D.from_FourierRealSeries(aux_nume_p, order)
    nume_m = Functions_1D.from_FourierRealSeries(aux_nume_m, order)
    deno = Functions_1D.from_FourierRealSeries(ww_sq, order)

    return nume_p * deno.inverse(), nume_m * deno.inverse()

    
    


# In[ ]:


def palipoly_coeffs_constructor(speed, coeffs):
    """
    Returns the palindromic coefficient vector of length 2*N+1 for
        P(z) = z^N [ speed + sum_{j=1}^{N} coeffs[j-1] * (z^j + z^{-j}) / 2],
    the Laurent polynomial representation of speed + w(x) at z = e^{ix}.
    """
    N = len(coeffs)
    polynom = zero_vector(RBF,2*N+1)
    polynom[N] = speed
    
    for jj, coef in enumerate(coeffs, start=1):
        polynom[N - jj] = coef * ONE_DIV_2
        polynom[N + jj] = coef * ONE_DIV_2

    return polynom


# In[ ]:


def big_bounds_compute_fast(polynom):
    """
    We want to compute upper bounds for
    dx^k P(z), with |z| less than 1.
    
    A first option would be using polynom_eval. Indeed:
        from explicit_funcs import polynom_eval
        big_bounds = vector(RBF, rk)
        for idx_der in range(rk):
            big_bounds[idx_der] = explicit_funcs.polynom_eval(ONE, abs_coeffs, idx_der=idx_der)
    
    However, this implementation is faster.
    """
    rk = len(polynom) - 1
    
    abs_coeffs = [coef.abs() for coef in polynom]
    
    big_bounds = vector(RBF, rk)
    for idx_der in range(rk):
        # For d>0, update weights w[k] *= (k-(idx_der-1)) for k>=idx_der
        bound = ZERO

        if idx_der == 0:
            bound = sum(abs_coeffs)

        else:
            mul = idx_der - 1
            for k in range(idx_der, rk):
                abs_coeffs[k] *= (k - mul)
                bound += abs_coeffs[k]

        big_bounds[idx_der] =  bound

    return big_bounds


# In[ ]:


def polynom_eval(x, coefficients, idx_der=0):
    """
    Evaluates P^{(idx_der)}(x) where P(x) = sum_k coefficients[k] * x^k.
    Uses the falling factorial for the derivative prefactor.
    """

    result = ZERO
    for idx_coef, coef in enumerate(coefficients):
        
        if idx_coef >= idx_der:
            #######################################
            # Compute k*(k-1)*...*(k-idx_der+1)
            aux_term = ONE
            
            for jj in range(idx_der):
                aux_term *= idx_coef-jj

            term = aux_term * coef * x**(idx_coef-idx_der)
            result += term

    return result


# In[ ]:


def root_polynomial_is_enclosed(root, polynom, rad, big_bounds):
    """
    Returns True if the ball B(root, rad) is a proven enclosure of a root of P.

    Uses the sufficient conditions:
        (|P(root)| + bound1) < rad * |P'(root)|_lower
        bound2 < |P'(root)|_lower
    where bound1, bound2 are remainder bounds built from big_bounds[k] >= sup|P^{(k)}|
    for k >= 3 over the ball.
    """

    #######################################
    # Remainder bounds built from derivatives >= 3
    bound1 = ZERO
    bound2 = ZERO

    for kk, val in enumerate(big_bounds,start=1):
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
    


# In[ ]:


def _vs_gen(speed, coeffs):
    """Returns vs = [speed, coeffs[0], ..., coeffs[N]] as an RBF vector."""
    vs = vector(RBF, len(coeffs) + 1)
    vs[0] = speed
    vs[1:] = coeffs
    return vs


# In[ ]:


def compute_resis_exis(speed, coeffs, ode_values):
    """
    Computes the existence residual vector of length n = len(vs),
    accumulating the ODE integral contributions over all intervals J in ode_values
    and adding the I6 norm term gk_exis_norm_sq.
    """
    vs = _vs_gen(speed, coeffs)

    n = len(vs)
    nm1 = n - 1
    poly_order = poly_ord(ode_values)

    D, S = precompute_DS(vs, vs)
    
    residues = vector(RBF, n)
    for counter, (J, values) in enumerate(ode_values, start=1):
        _, _, x0, rr = interval_from_J(J)
        pre = pre_resis_exis(x0, rr, vs, D, S, poly_order)

        for idx, func_values in enumerate(values):
            kk = idx + 1
            residues[idx] += I_sols_resi_exis_from_precomp(func_values, kk, pre, nm1)

        print_iter_J(counter, text=f"{J}")

    residues += gk_exis_norm_sq(vs)
    return residues
    

def compute_resis_stab(lamb, the, speed, coeffs, ode_values):
    """
    Computes the stability residual vector of length 4*(n-1)+2,
    accumulating ODE integral contributions over all intervals J
    and adding the I6 norm term gk_stab_norm_sq.
    """
    vs = _vs_gen(speed, coeffs)
    
    nm1 = len(vs) - 1
    nn = 4*nm1 + 2
    poly_order = poly_ord(ode_values)
    
    Wp0, Wm0, DS_data = DS_stab_helper(vs, lamb, the)
    
    residues = vector(RBF, nn)
    for counter, (J, values) in enumerate(ode_values, start=1):        
        _, _, x0, rr = interval_from_J(J)
        pre = pre_resis_stab(x0, rr, vs, Wp0, Wm0, DS_data, the, poly_order)
    
        for idx, func_values in enumerate(values):
            kk = idx + 1
            residues[idx] += I_sols_resi_stab_from_precomp(func_values, kk, pre, nm1)

        print_iter_J(counter, text=f"{J}")
    
    residues += gk_stab_norm_sq(vs, the)
    return residues


def compute_resis_Jfv(lamb, the, fv, speed, coeffs, ode_values):
    """
    Computes the J_{f,v} residual pair (resi_p, resi_m),
    accumulating contributions from ode_values and adding the L2 norm of fv split as (fvp, fvm).
    """
    vs = _vs_gen(speed, coeffs)

    poly_order = poly_ord(ode_values)
    nm1 = len(fv) // 2
    fvp, fvm = fv[:nm1], fv[nm1:]
    
    Wp0, Wm0, DS_data = DS_stab_helper(vs, lamb, the)

    resi_p = ZERO; resi_m = ZERO
    for counter, (J, values) in enumerate(ode_values, start=1):
        _, _, x0, rr = interval_from_J(J)
        pre = pre_resis_stab_Jfv(x0, rr, vs, Wp0, Wm0, DS_data, the, poly_order, fvp, fvm)

        resi_p += I_sols_resi_Jfv_from_precomp(values[0], 0, pre, nm1)
        resi_m += I_sols_resi_Jfv_from_precomp(values[1], 1, pre, nm1)
        
        print_iter_J(counter, text=f"{J}")

    resi_p += auxiliar_funcs.norm_vector_sq(fvp) * TWOPI
    resi_m += auxiliar_funcs.norm_vector_sq(fvm) * TWOPI
    return resi_p, resi_m


# In[ ]:


def _build_CSvals(x0, rr, n, rk):
    """
    Builds Cvals[ll][idx_jj] and Svals[ll][idx_jj] for jj in [-N, 2N]
    and ll in [0, 2rk], where idx_jj is the position in that range.
    """
    max_Nk = 2*rk
    max_Nj = 2*n-1
    min_Nj = n-1
    
    Cvals = [ [] for _ in range(max_Nk + 1) ] # Maybe is a good idea change this from list to dict
    Svals = [ [] for _ in range(max_Nk + 1) ]
    
    jjs_cs   = list(range(-min_Nj, max_Nj))
    for jj in jjs_cs:
        C, S = auxiliar_funcs.build_CS_for_jj(jj, max_Nk, x0=x0, rr=rr)
        for ll in range(max_Nk + 1):
            Cvals[ll].append(C[ll])
            Svals[ll].append(S[ll])

    return Cvals, Svals

def _build_Evals(x0, rr, n, rk):
    """
    Builds Evals[ll] as a dict {jj: C[ll] - i*S[ll]} for jj in [-(2n+1), 2n+1]
    and ll in [0, 2rk]. Used for the stability integrals II4, II5.
    """
    max_Nk = 2*rk
    max_Nj = 2*n + 2
    min_Nj = 2*n + 1
    
    Evals = [ {} for _ in range(max_Nk + 1) ]

    jjs_eval = list(range(-min_Nj, max_Nj))   # for Evals dict
    for jj in jjs_eval:
        C, S = auxiliar_funcs.build_CS_for_jj(jj, max_Nk, x0=x0, rr=rr)
        for ll in range(max_Nk + 1):
            Evals[ll][jj] = C[ll] - I * S[ll]

    return Evals
    


# In[ ]:


def pre_resis_exis(x0, rr, vs, Dvs, Svs, rk):
    """Precomputes the II1..II5 integral tables for the existence residual on interval [x0-rr, x0+rr]."""
    n = len(vs)
    rk1 = rk + 1
    
    Cvals, Svals = _build_CSvals(x0, rr, n, rk)
    
    II1 = [ [II1_exis(l1, l2,     n, Cvals          ) for l2 in range(rk1)] for l1 in range(rk1)]
    II2 = [ [II2_exis(l1, l2,     n, Cvals, Dvs, Svs) for l2 in range(rk1)] for l1 in range(rk1)]
    II3 = [ [II3_exis(l1, l2, vs, n, Cvals          ) for l2 in range(rk1)] for l1 in range(rk1)]

    II4 = [II4_exis(ll, vs, n, Svals, Cvals) for ll in range(rk1)]
    II5 = [II5_exis(ll, vs, n, Svals, Cvals) for ll in range(rk1)]

    pack = {"II1": II1, "II2": II2, 
            "II3": II3, "II4": II4, 
            "II5": II5}
    return pack


def pre_resis_stab(x0, rr, vs, Wp0, Wm0, DS_data, the, rk):
    """Precomputes the II1..II5 integral tables for the stability residual on interval [x0-rr, x0+rr]."""
    n = len(vs)
    rk1 = rk + 1
    
    Cvals, Svals = _build_CSvals(x0, rr, n, rk)
    Evals = _build_Evals(x0, rr, n, rk)

    II1 = [[II1_stab(l1, l2, n, Cvals, DS_data) for l2 in range(rk1)] for l1 in range(rk1)]
    II2 = [[II2_stab(l1, l2, n, Cvals, DS_data) for l2 in range(rk1)] for l1 in range(rk1)]
    II3 = [[II3_stab(l1, l2, n, Cvals, DS_data) for l2 in range(rk1)] for l1 in range(rk1)]

    II4 = [II4_stab(ll, vs, Wp0, Wm0, the, Evals) for ll in range(rk1)]
    II5 = [II5_stab(ll, vs, the, Evals) for ll in range(rk1)]

    pack = {"II1": II1, "II2": II2, 
            "II3": II3, "II4": II4, 
            "II5": II5}
    return pack


def pre_resis_stab_Jfv(x0, rr, vs, Wp0, Wm0, DS_data, the, rk, fvp, fvm):
    """Precomputes the II1..II3 and II45 integral tables for the J_{f,v} residual on [x0-rr, x0+rr]."""
    n = len(vs)
    rk1 = rk + 1
    
    Cvals, Svals = _build_CSvals(x0, rr, n, rk)
    Evals = _build_Evals(x0, rr, n, rk)

    II1 = [ [II1_stab(l1, l2, n, Cvals, DS_data) for l2 in range(rk1)] for l1 in range(rk1)]
    II2 = [ [II2_stab(l1, l2, n, Cvals, DS_data) for l2 in range(rk1)] for l1 in range(rk1)]
    II3 = [ [II3_stab(l1, l2, n, Cvals, DS_data) for l2 in range(rk1)] for l1 in range(rk1)]

    II45 = I45_stab_fv_pm(Evals, vs, Wp0, Wm0, the, fvp, fvm)

    pack = {"II1": II1, "II2": II2,
            "II3": II3, "II45": II45}
    return pack


# In[ ]:


def precompute_DS(w1, w2):
    """
    Computes the difference and sum convolution arrays of w1 and w2:
        D[j1-j2+(n-1)] += w1[j1] * conj(w2[j2])
        S[j1+j2]       += w1[j1] * conj(w2[j2])
    Both D and S have length 2n-1. w2 is cast to CBF internally.
    """
    if len(w1) != len(w2):
        raise ValueError("`w1` and `w2` must have the same length.")
        
    n = len(w1)
    D = [ZERO] * (2*n-1)   # differences
    S = [ZERO] * (2*n-1)   # sums

    for j1, v1 in enumerate(w1):
        for j2, v2 in enumerate(w2):
            prod = v1 * CBF(v2).conjugate()

            D[j1 - j2 + (n-1)] += prod
            S[j1 + j2] += prod
            
    return D, S

def DS_stab_helper(vs, la, the):
    """
    Builds Wp, Wm and precomputes all DS pairs needed for the stability integrals.
        Wp[0] = i + la - i*the*vs[0],   Wp[k] = -i*the*vs[k]  for k >= 1
        Wm[0] = i - la - i*eht*vs[0],   Wm[k] = -i*eht*vs[k]  for k >= 1
    Returns (Wp0, Wm0, DS) where DS is a dict keyed by pairs like ("Wp","vs").
    """
    eht = 1 - the

    Wp = -I * the * vs
    Wm = -I * eht * vs

    Wp0 = I + la - I * the * vs[0]
    Wm0 = I - la - I * eht * vs[0]
    Wp[0] = Wp0; Wm[0] = Wm0

    vecs = {"vs": vs, "Wp": Wp, "Wm": Wm}

    pairs = [
        ("vs", "vs"),
        ("Wp", "vs"),
        ("Wm", "vs"),
        ("Wp", "Wp"),
        ("Wm", "Wm"),
    ]

    DS = {}
    for a, b in pairs:
        D, S = precompute_DS(vecs[a], vecs[b])
        DS[(a, b)] = (D, S)

    return Wp0, Wm0, DS


# In[ ]:


def gk_exis_norm_sq(vs):
    """
    Computes the existence I6 norm vector: I6[k-1] = pi * k^2 * sum_{m>=k} vs[m]^2 (with vs[k]^2/2 for m=k).
    The last entry is always zero (asserted).
    """
    I6 = vector(RBF, len(vs))
    for kk, vk in enumerate(vs):
        if kk==0:
            continue
        term = vk**2 * ONE_DIV_2
        
        for mm in range(1, len(vs)-kk):
            term += vs[mm+kk]**2
            
        I6[kk-1] = term*kk**2

    assert I6[-1] == ZERO
    return I6*PI

def gk_stab_norm_sq(vs, the):
    """
    Computes the stability I6 norm vector of length 4*(n-1)+2, with blocks:
        [I6_the | I6_eht | I6_eht | I6_the | 0 | 0]
    where I6_the[k-1] = (pi/2) * (k-the)^2 * sum_{m>=k} vs[m]^2.
    The last two entries are always zero (asserted).
    """

    eht = 1 - the
    
    n = len(vs)
    nm1 = n - 1
    nn = 4*nm1 + 2

    I6 = vector(RBF, nn)
    
    I6_the = vector(RBF, nm1)
    I6_eht = vector(RBF, nm1)
    for kk in range(1, n):
        kkm1 = kk - 1
        
        term = ZERO
        for mm in range(n - kk):
            term += vs[mm + kk]**2
            
        I6_the[kkm1] = term * (kk - the)**2
        I6_eht[kkm1] = term * (kk - eht)**2

    coef = PI * ONE_DIV_2
    I6[      :   nm1] = I6_the[:]
    I6[  nm1 : 2*nm1] = I6_eht[:]
    I6[2*nm1 : 3*nm1] = I6_eht[:]
    I6[3*nm1 : 4*nm1] = I6_the[:]

    assert I6[-1] == ZERO
    assert I6[-2] == ZERO
    return I6 * coef


# In[ ]:


def I_sols_resi_stab_from_precomp(values, kk, pre, nm1):
    """Assembles the kk-th stability residual component from precomputed II1..II5 tables."""
    II1, II2 = pre["II1"], pre["II2"]
    II3, II4 = pre["II3"], pre["II4"]
    II5      = pre["II5"]
        
    kkm1 = kk - 1

    sJ = 0 if (kk <= 2*nm1 or kk == 4*nm1 + 1) else 1    
    
    residue = ZERO
    for l1, alpha1 in enumerate(values):
        
        if kk <= 4*nm1:
            new_I4 = (alpha1 * CBF(II4[l1][kkm1])).real()
            new_I5 = (alpha1 * CBF(II5[l1][kkm1])).real()
      
            residue += new_I4 + new_I5
          
        
        for l2, alpha2 in enumerate(values):
            alpha1_alpha2 = alpha1 * alpha2.conjugate()
            
            new_I1 =  alpha1_alpha2 * II1[l1][l2][sJ]
            new_I2 =  alpha1_alpha2 * II2[l1][l2]        
            new_I3 = (alpha1_alpha2 * II3[l1][l2][sJ]).real()
                          
            residue += new_I1 + new_I2 + new_I3
    
    return residue.real()


def I_sols_resi_exis_from_precomp(values, kk, pre, nm1):
    """Assembles the kk-th existence residual component from precomputed II1..II5 tables."""
    II1, II2 = pre["II1"], pre["II2"]
    II3, II4 = pre["II3"], pre["II4"]
    II5      = pre["II5"]
        
    kkm1 = kk - 1

    residue = ZERO
    for l1, alpha1 in enumerate(values):
        
        if kk <= nm1:
            new_I4 = (alpha1 * II4[l1][kkm1]).real()
            new_I5 = (alpha1 * II5[l1][kkm1]).real()
      
            residue += new_I4 + new_I5
        
        for l2, alpha2 in enumerate(values):
            alpha1_alpha2 = alpha1 * alpha2.conjugate()
            
            new_I1 =  alpha1_alpha2 * II1[l1][l2]
            new_I2 =  alpha1_alpha2 * II2[l1][l2]        
            new_I3 = (alpha1_alpha2 * II3[l1][l2]).real()
                          
            residue += new_I1 + new_I2 + new_I3

    return residue.real()


def I_sols_resi_Jfv_from_precomp(values, sJ, pre, nm1):
    """Assembles the sJ-th J_{f,v} residual component from precomputed II1..II3 and II45 tables."""
    II1, II2  = pre["II1"], pre["II2"]
    II3, II45 = pre["II3"], pre["II45"]
    
    residue = ZERO
    for l1, alpha1 in enumerate(values):
        
        new_I45 = (alpha1 * CBF(II45[l1][sJ])).real()
        
        residue += new_I45
        
        for l2, alpha2 in enumerate(values):
            alpha1_alpha2 = alpha1 * alpha2.conjugate()
            
            new_I1 =  alpha1_alpha2 * II1[l1][l2][sJ]
            new_I2 =  alpha1_alpha2 * II2[l1][l2]        
            new_I3 = (alpha1_alpha2 * II3[l1][l2][sJ]).real()
            
            residue += new_I1 + new_I2 + new_I3
    
    return residue.real()


# In[ ]:


def II1_exis(l1, l2, n, Cvals):
    p = l1 + l2
    shift = n - 1
    return Cvals[p][shift]

def II2_exis(l1, l2, n, Cvals, D, S):
    if l1 == 0 or l2 == 0:
        return ZERO

    Ck = Cvals[l1 + l2 - 2]
    shift = n - 1

    acc = ZERO
    for tt in range(2*n - 1):
        acc += D[tt] * Ck[tt] + S[tt] * Ck[tt + shift]

    return acc * (l1 * l2 * ONE_DIV_2)

def II3_exis(l1, l2, vs, n, Cvals):
    if l1 == 0:
        return ZERO

    Ck = Cvals[l1 + l2 - 1]
    shift = n - 1

    acc = ZERO
    for jj, vj in enumerate(vs):
        acc += vj * Ck[jj + shift]

    return 2 * I * l1 * acc

def II4_exis(ll, vs, n, Svals, Cvals):
    shift = n - 1
    Srow = Svals[ll]
    C0 = Cvals[ll][shift]   

    P = [ZERO] * n
    for t in range(2, n):    
        vt = vs[t]
        base = shift + t
        for k in range(1, t):
            # m = t-k >= 1
            P[k] += vt * Srow[base - k]  

    row = [ZERO] * (n - 1)
    for k in range(1, n):
        row[k - 1] = (k * vs[k] * C0) + (2 * k * I) * P[k]

    return row


def II5_exis(l, vs, n, Svals, Cvals):
    row = [ZERO] * (n - 1)
    if l < 1:
        return row

    shift = n - 1

    # --- T = sum_j twap_j * C_J^{l-1,j}
    Cprev = Cvals[l - 1]
    T = ZERO
    for j in range(n):
        T += vs[j] * Cprev[shift + j]

    # --- Precompute III_m for m=1..N 
    Sp = Svals[l - 1]
    III_m = [ZERO] * n
    for m in range(1, n):
        III_m[m] = III5_exis(l, m, vs, n, Sp)

    # --- For each k, compute:
    # term1 = i*l*k*twap_k*T
    # term2 = -l*k * sum_{m=1}^{N-k} twap_{m+k} * III_m[m]
    for k in range(1, n):
        idx = k - 1

        acc = ZERO
        # m = 1..(n-1-k)
        for m in range(1, n - k):
            acc += vs[m + k] * III_m[m]

        row[idx] = (I * l * k) * vs[k] * T - (l * k) * acc

    return row


def III5_exis(l, m, vs, n, Sprev):
    shift = n - 1
    base = shift + m

    acc = ZERO
    for jj, vj in enumerate(vs):
        acc += vj * (Sprev[base + jj] + Sprev[base - jj])
    return acc


# In[ ]:


def II1_stab(l1, l2, n, Cvals, DS_data):
    (Dp, Sp) = DS_data[("Wp", "Wp")]
    (Dm, Sm) = DS_data[("Wm", "Wm")]
    
    p = l1 + l2
    Ck = Cvals[p]
    shift = n - 1

    acc_p = ZERO
    acc_m = ZERO

    for tt in range(2*n - 1):
        acc_p += Dp[tt] * Ck[tt] + Sp[tt] * Ck[tt + shift]
        acc_m += Dm[tt] * Ck[tt] + Sm[tt] * Ck[tt + shift]

    return [ONE_DIV_2 * acc_p, ONE_DIV_2 * acc_m]


def II2_stab(l1, l2, n, Cvals, DS_data):
    (D, S) = DS_data[("vs", "vs")]
    return II2_exis(l1, l2, n, Cvals, D, S)


def II3_stab(l1, l2, n, Cvals, DS_data):
    if l2 == 0:
        return [ZERO]*2

    (Dp, Sp) = DS_data[("Wp", "vs")]
    (Dm, Sm) = DS_data[("Wm", "vs")]

    p = l1 + l2 - 1
    Ck = Cvals[p]
    shift = n - 1

    acc_p = ZERO
    acc_m = ZERO
    for tt in range(2*n - 1):
        acc_p += Dp[tt] * Ck[tt] + Sp[tt] * Ck[tt + shift]
        acc_m += Dm[tt] * Ck[tt] + Sm[tt] * Ck[tt + shift]

    coef = - l2
    return [coef * acc_p, coef * acc_m]



def II4_stab(ll, vs, Wp0, Wm0, the, Evals):

    eht = 1 - the

    Erow = Evals[ll]      # E^{l,...}
    n = len(vs)
    coef = -I * ONE_DIV_2

    alpha_p = -I * the
    alpha_m = -I * eht

    t0_sp, t0_sm, T_sp, T_sm = III4_stab(vs, Erow)

    res_pp = [ZERO] * (n - 1)
    res_pm = [ZERO] * (n - 1)
    res_mp = [ZERO] * (n - 1)
    res_mm = [ZERO] * (n - 1)

    for kk in range(1, n):
        idx = kk - 1
        mmax = n - kk

        # dot products against w_m = vs[m+kk]
        dot_t0_sp = ZERO
        dot_T_sp  = ZERO
        dot_t0_sm = ZERO
        dot_T_sm  = ZERO

        for m in range(mmax):
            w = vs[m + kk]
            dot_t0_sp += w * t0_sp[m]
            dot_T_sp  += w * T_sp[m]
            dot_t0_sm += w * t0_sm[m]
            dot_T_sm  += w * T_sm[m]

        acc_pp = Wp0 * dot_t0_sp + alpha_p * dot_T_sp
        acc_pm = Wp0 * dot_t0_sm + alpha_p * dot_T_sm
        acc_mp = Wm0 * dot_t0_sp + alpha_m * dot_T_sp
        acc_mm = Wm0 * dot_t0_sm + alpha_m * dot_T_sm

        kt = kk - the
        ke = kk - eht

        res_pp[idx] = coef * acc_pp * kt
        res_pm[idx] = coef * acc_pm * ke
        res_mp[idx] = coef * acc_mp * ke
        res_mm[idx] = coef * acc_mm * kt

    return res_pp + res_pm + res_mp + res_mm


def III4_stab(vs, Erow):
    n = len(vs)

    t0_sp = [ZERO] * n
    t0_sm = [ZERO] * n
    T_sp  = [ZERO] * n
    T_sm  = [ZERO] * n

    for m in range(n-1):
        sp =  m
        sm = -(m + 1)

        t0_sp[m] = 2 * Erow[sp]
        t0_sm[m] = 2 * Erow[sm]

        acc_sp = ZERO
        acc_sm = ZERO
        for jj, vj in enumerate(vs):
            if jj == 0:
                continue
            acc_sp += vj * (Erow[sp - jj] + Erow[sp + jj])
            acc_sm += vj * (Erow[sm - jj] + Erow[sm + jj])

        T_sp[m] = acc_sp
        T_sm[m] = acc_sm

    return t0_sp, t0_sm, T_sp, T_sm

    
def II5_stab(ll, vs, the, Evals):

    n = len(vs)
    nm1 = n - 1

    if ll == 0:
        return [ZERO] * (4*nm1)

    eht = 1 - the
    Erow = Evals[ll - 1]

    III_p = [ZERO] * n
    III_m = [ZERO] * n

    for mm in range(n):
        shift_p =  mm
        shift_m = -(mm + 1)

        III_p[mm] = III5_stab(vs, shift_p, Erow)
        III_m[mm] = III5_stab(vs, shift_m, Erow)

    res_pp = [ZERO] * nm1
    res_pm = [ZERO] * nm1
    res_mp = [ZERO] * nm1
    res_mm = [ZERO] * nm1

    coef = I * ONE_DIV_2 * ll

    for kk in range(1, n):
        idx = kk - 1
        
        acc_p = ZERO
        acc_m = ZERO

        for mm in range(n-kk):
            w = vs[mm + kk]
            acc_p += w * III_p[mm]
            acc_m += w * III_m[mm]

        kt = (kk - the)  
        ke = (kk - eht)  

        res_pp[idx] = coef * acc_p * kt
        res_pm[idx] = coef * acc_m * ke
        res_mp[idx] = coef * acc_p * ke
        res_mm[idx] = coef * acc_m * kt

    return res_pp + res_pm + res_mp + res_mm


def III5_stab(vs, shift, Erow_lm1):
    s = shift
    
    acc = ZERO
    for jj, vj in enumerate(vs):
        acc += vj * (Erow_lm1[s - jj] + Erow_lm1[s + jj])
    return acc


# In[1]:


def I45_stab_fv_pm_l(E_l, E_prev, l, vs, Wp0, Wm0, the, fvp, fvm):
    n = len(vs)
    nm1 = n - 1

    alpha_p = -I * the
    alpha_m = -I * (1 - the)

    def build_conv(E, vs, W0, alpha):
        W = vector(CBF, [W0] + [alpha * term for term in vs[1:]])
        
        out = vector(CBF, nm1)
        for m in range(nm1):
            acc = ZERO
            for j, wj in enumerate(W):
                acc += wj * (E[m - j] + E[m + j])
            out[m] = acc
        return out

    # ----- II4 -----
    CWp = build_conv(E_l, vs, Wp0, alpha_p)
    CWm = build_conv(E_l, vs, Wm0, alpha_m)

    II4_p = ZERO
    for (fm, cw) in zip(fvp, CWp):
        II4_p += fm.conjugate() * cw

    II4_m = ZERO
    for (fm, cw) in zip(fvm, CWm):
        II4_m += fm.conjugate() * cw

    # ----- II5 -----
    if l == 0:
        II5_p = ZERO
        II5_m = ZERO
    else:

        III = vector(CBF, nm1)
        for m in range(nm1):
            acc = ZERO
            for j, vj in enumerate(vs):
                acc += vj * (E_prev[m - j] + E_prev[m + j])
            III[m] = acc

        II5_p = ZERO
        for (fm, Im) in zip(fvp, III):
            II5_p += fm.conjugate() * Im
        II5_p *= -l

        II5_m = ZERO
        for (fm, Im) in zip(fvm, III):
            II5_m += fm.conjugate() * Im
        II5_m *= -l

    return [II4_p + II5_p, II4_m + II5_m]


def I45_stab_fv_pm(Evals, vs, Wp0, Wm0, the, fvp, fvm):
    out = []
    max_Nk_for_E = (len(Evals) + 1) // 2
    for l, E_l in enumerate(Evals[:max_Nk_for_E]):
        E_prev = None if l == 0 else Evals[l - 1]
        out.append(I45_stab_fv_pm_l(E_l, E_prev, l, vs, Wp0, Wm0, the, fvp, fvm))
    return out



# In[ ]:


def construct_exis_mat(ode_values):
    """
    Assembles the existence matrix (n x n, CBF) by accumulating E*F over all intervals J,
    dividing by 2*pi, and subtracting 1 on the first subdiagonal.
    """
    n = len(ode_values[0][1]) # Maybe is a good idea construct this as a new class. n, poly_order would be attributes.
    poly_ord_p1 = poly_ord(ode_values) + 1

    exis_mat = matrix(CBF, n)   # (n x n) zero matrix
    for counter, (J, values) in enumerate(ode_values, start=1):
        expo_vals = auxiliar_funcs.precompute_expo_vals(J, 0, n, poly_ord_p1)
    
        F = matrix(CBF, values).transpose() # rows = ll, cols = kk
        E = matrix(CBF, expo_vals)          # rows = jj, cols = ll
        exis_mat += E * F                   # rows = jj, cols = kk
        
        print_iter_J(counter, text=f"{J}")

    exis_mat /= TWOPI    
    for kk in range(n-1):
        exis_mat[kk+1, kk] -= ONE

    return exis_mat



# In[ ]:


def construct_stab_mat(ode_values):
    """
    Assembles the stability matrix (2n x 2n, CBF) by accumulating block-permuted E*F
    over all intervals J, dividing by 2*pi, and `subtracting 1 on two subdiagonals`.
    """

    nn = len(ode_values[0][1])
    N = (nn-2) // 4
    n = N + 1
    poly_ord_p1 = poly_ord(ode_values) + 1

    stab_mat = matrix(CBF, 2*n)   # 2n square matrix
    for counter, (J, values) in enumerate(ode_values, start=1):
        expo_vals = auxiliar_funcs.precompute_expo_vals(J, -1, N, poly_ord_p1)

        F = matrix(CBF, values).transpose() # rows = ll, cols = kk
        E = matrix(CBF, expo_vals)          # rows = jj, cols = ll
        proj_mat = E * F         

        # Top block rows 0..N, cols 0..2N (last col of stab_mat is Pi^-)
        stab_mat[0:n, 0:N]   += proj_mat[:, N:2*N]
        stab_mat[0:n, N:2*N] += proj_mat[:, 0:N]    
        stab_mat[0:n, 2*N]   += proj_mat[:, 4*N]

        # Bottom block rows N+1..2N+1
        stab_mat[n:2*n, 0:2*N] += proj_mat[:, 2*N:4*N]    
        stab_mat[n:2*n, 2*N+1] += proj_mat[:, 4*N+1]   # mhom (Pi^-)

        print_iter_J(counter, text=f"{J}")

    stab_mat /= TWOPI
    for kk in range(N):
        stab_mat[  kk+1,   kk] -= ONE
        stab_mat[n+kk+1, N+kk] -= ONE

    return stab_mat


# In[ ]:


def bfv_constructor(ode_values, nm1, ord_mat):
    poly_ord_p1 = poly_ord(ode_values) + 1
    n = nm1 + 1
    
    bfv = vector(CBF, ord_mat)
    for counter, (J, values) in enumerate(ode_values, start=1):
        expo_vals = auxiliar_funcs.precompute_expo_vals(J, -1, nm1, poly_ord_p1)

        F = matrix(CBF, values).transpose() # rows = ll, cols = 4
        E = matrix(CBF, expo_vals)          # rows = jj, cols = ll
        proj_mat = E * F
            
        bfv[:n] += proj_mat.column(0)
        bfv[n:] -= proj_mat.column(1) # <<<<-------- THE SIGN IS IMPORTANT

        print_iter_J(counter, text=f"{J}")
        
    return bfv / TWOPI


# In[ ]:


def mat_exis_add_radii(M, bounds):
    """
    Returns the Frobenius norm of the Weyl perturbation bound matrix with entries
        result[i,j] = r_i*||col_j|| + r_j*||col_i|| + r_i*r_j,
    where r_i = (bounds[i] / 2*pi)^{1/2}.
    """
    ord_mat = len(bounds)
    
    #######################################
    # Basic dimension check: M must be square
    verify.square_mat_order(M, ord_mat)

    #######################################
    # Precompute auxiliary quantities        
    res_terms = div_twopi_sqrt_vec(bounds) # bounds contains squared L2 integrals, we normalize them and remove the square.
    mat_terms = col_norms(M)
    
    #######################################
    # Fill the bound matrix
    result = matrix(CBF, ord_mat)
    for ii in range(ord_mat):
        for jj in range(ord_mat):
            term1 = res_terms[ii] * mat_terms[jj]
            term2 = res_terms[jj] * mat_terms[ii]
            term3 = res_terms[ii] * res_terms[jj]

            result[ii, jj] = term1 + term2 + term3
    
    #######################################
    # Return Frobenius norm bound (scalar)
    return auxiliar_funcs.fro_norm(result)


# In[ ]:


def mat_stab_add_radii(stab_mat, ode_bounds):
    """
    Same as mat_exis_add_radii but for the stability matrix, where the residuals
    are not one-to-one with columns. Reassigns ode_bounds into regular_resii
    following the block structure of stab_mat before calling mat_exis_add_radii.

    Here the residuals are not one-to-one with the columns.
    Denoting by b := ode_bounds, and bk_{j} the j-th Fourier coefficient of bk.
    The structure with len(coeffs) = 2 would be the following: ord_mat = 10 and 
        |b3_{-1} b4_{-1} b1_{-1} b2_{-1} b9_{-1} 0       |
        |b3_{ 0} b4_{ 0} b1_{ 0} b2_{ 0} b9_{ 0} 0       |
        |b3_{ 1} b4_{ 1} b1_{ 1} b2_{ 1} b9_{ 1} 0       |
        |b5_{-1} b6_{-1} b7_{-1} b8_{-1} 0       b10_{-1}|
        |b5_{ 0} b6_{ 0} b7_{ 0} b8_{ 0} 0       b10_{ 0}|
        |b5_{ 1} b6_{ 1} b7_{ 1} b8_{ 1} 0       b10_{ 1}|
         r1      r2      r3      r4      r5      r5
    We call `mat_exis_add_radii()` with the corrected residuals: `regular_resii` using:
        rk^2 = sum_{l=-1}^1 |bk1_l|^2 + |bk2_l|^2 <= sum_{l=-1}^1 (bk1 + bk2) / (2 pi), where bk1 = ode_bounds(k1).
    Normalization and removing the square is done in mat_exis_add_radii.
    """

    # Compute `regular_resii` and call `mat_exis_add_radii`
    regular_resii = aux_sort_stab_mat(ode_bounds)
    return mat_exis_add_radii(stab_mat, regular_resii)


# In[ ]:


def compute_hk_norms_exis_ap(ode_values):
    """Computes the different quantities required in ... with the provided approximations""" 
    n = len(ode_values[0][1])
    N = n - 1
    poly_ord_p1 =  poly_ord(ode_values) + 1

    hk_norm_ap_0 = vector(RBF, N)
    hk_norm_ap_1 = vector(RBF, N)

    for (J, values) in ode_values:
        _, _, _, rr = interval_from_J(J)
        terms = _precomp_hk_norms_ap(rr, 2*poly_ord_p1-1)

        for ik, vk in enumerate(values[:-1]): # Pi already covered so we only return the hk, for k in [1,N]
            hk_norm_ap_0[ik] += _hk_norms_0(vk, terms)
            hk_norm_ap_1[ik] += _hk_norms_1(vk, terms)

    return hk_norm_ap_0, hk_norm_ap_1


# In[ ]:


def compute_hk_norms_stab_ap(ode_values):
    """Computes the different quantities required in ... with the provided approximations""" 
    nn = len(ode_values[0][1])
    N = (nn-2) // 4
    n = N + 1
    poly_ord_p1 =  poly_ord(ode_values) + 1

    hk_norm_ap_p0 = vector(RBF, 2*N+1)
    hk_norm_ap_m0 = vector(RBF, 2*N+1)
    hk_norm_ap_p1 = vector(RBF, 2*N+1)
    hk_norm_ap_m1 = vector(RBF, 2*N+1)

    for (J, values) in ode_values:
        _, _, _, rr = interval_from_J(J)
        terms = _precomp_hk_norms_ap(rr, 2*poly_ord_p1-1)

        shift_pm = 2*N
        for ik in range(shift_pm):
            vals_kp = values[ik]
            vals_km = values[ik + shift_pm]
            hk_norm_ap_p0[ik] += _hk_norms_0(vals_kp, terms)
            hk_norm_ap_m0[ik] += _hk_norms_0(vals_km, terms)
            hk_norm_ap_p1[ik] += _hk_norms_p(vals_kp, terms)
            hk_norm_ap_m1[ik] += _hk_norms_m(vals_km, terms)

        last_vals_p = values[-2]
        last_vals_m = values[-1]
        hk_norm_ap_p0[-1] += _hk_norms_0(last_vals_p, terms)
        hk_norm_ap_m0[-1] += _hk_norms_0(last_vals_m, terms)
        hk_norm_ap_p1[-1] += _hk_norms_p(last_vals_p, terms)
        hk_norm_ap_m1[-1] += _hk_norms_m(last_vals_m, terms)

    return hk_norm_ap_p0, hk_norm_ap_m0, hk_norm_ap_p1, hk_norm_ap_m1


# In[ ]:


# Some secondary functions, not enough general to be placed in auxiliar_funcs
def div_twopi_sqrt_vec(v):
    """Only RBF"""
    verify.oo(v, len(v))
    half = ONE_DIV_2
    twopi = TWOPI
    return vector(RBF, [(x / twopi)**half for x in v])


def col_norms(M):
    """Returns a vector, k-th entry is the norm of the k-th column"""
    return vector(RBF, [col.norm() for col in M.columns()])
    

def aux_sort_stab_mat(ode_bounds):
    """Auxiliar sortering function, see mat_stab_add_radii docstring """
    nn = len(ode_bounds)
    nm1 = nn // 4
    ord_mat = 2*nm1 + 2

    regular_resii = vector(RBF, ord_mat)
    for kk in range(nm1):
        regular_resii[kk      ] = ode_bounds[kk + nm1] + ode_bounds[kk + 2*nm1]
        regular_resii[kk + nm1] = ode_bounds[kk      ] + ode_bounds[kk + 3*nm1]
        
    regular_resii[-2] = ode_bounds[-2]
    regular_resii[-1] = ode_bounds[-1]

    return regular_resii


def poly_ord(values):
    """Helper function to access the order of the approximating polynomials"""
    return len(values[0][1][0]) - 1


def _precomp_hk_norms_ap(rr, max_ord):
        terms = vector(RBF, max_ord)
        for ll in range(max_ord):
            if ll % 2 == 1:
                continue
            else:
                ll1 = ll + 1
                terms[ll] = 2 * rr**ll1 / RBF(ll1)
        return terms


def _hk_norms_0(vals, terms):
    acc = ZERO
    for ij, vj in enumerate(vals):
        for ik, vk in enumerate(vals):
            acc += vj * vk.conjugate() * terms[ij + ik]
    return acc.real()


def _hk_norms_1(vals, terms):
    vals_der = _vals_der(vals)
    return _hk_norms_0(vals_der, terms)


def _hk_norms_p(vals, terms):
    return _hk_norms_1(vals, terms)


def _hk_norms_m(vals, terms):
    vals_der = _vals_der(vals)
    vals_mod = I*vals + vals_der
    return _hk_norms_0(vals_mod, terms)


def _vals_der(vals):
    return vector(CBF, [idx * val for idx, val in enumerate(vals)][1:] + [ZERO])


# In[ ]:




