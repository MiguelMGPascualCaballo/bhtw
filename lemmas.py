#!/usr/bin/env python
# coding: utf-8

# In[10]:


from sage.all import *
from classes import FourierRealSeries, Functions_1D
import auxiliar_funcs
from methods import image_1D, integral_1D
from parameters import (
    RBF, CBF, 
    SQRT2, ONE_DIV_2, ZERO, ONE, TWO, PI, TWOPI,
    VERBOSE,
)
import load_data
import verify
from printing_macros import print_iter_N, print_abcde

BOUNDS, BOUNDS_STR = load_data.load_bounds()


# In[ ]:


def norm_xi_twap_sq(): # (Lemma 2.3)
    """
    Compute and validate an enclosure for the L2 norm of the residual xi
    associated with the approximate traveling wave.

    The residual is defined as

        xi = c * v' + Hv + v' * v,

    where c is the wave speed and v is the approximate traveling wave. The
    function evaluates xi symbolically, computes its norm in the Fourier-based
    convention used in the paper, and checks whether the resulting enclosure
    is below the bound prescribed in `BOUNDS['resi_exis_L2']`.

    The norm convention in the paper is expressed using only positive
    Fourier modes for real odd functions. If

        xi(x) = sum_{k in Z} xi_k e^{ikx},

    then the norm used in the paper is

        (sum_{k >= 1} |xi_k|^2)^(1/2).

    Since `xi.norm_homogeneous()` returns

        (sum_{k in Z-{0}} |xi_k|^2)^(1/2),

    and xi is assumed to be real and odd, we divide by `sqrt(2)` to
    recover the norm used in the argument.
    """

    #######################################
    # Output information
    lemma_label = "Lemma 2.3"
    bound_labels = ['resi_exis_L2']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]

    #######################################
    # Load data
    from explicit_funcs import equ_tw_eval_symbolic
    
    speed  = load_data.speed()
    coeffs = load_data.coeffs()
    
    norm_xiL2_save = BOUNDS['resi_exis_L2']
    
    #######################################
    # Compute xi
    xi = equ_tw_eval_symbolic(speed, coeffs)

    # Compute the norm
    norm_xiL2_comp = xi.norm_homogeneous() / SQRT2

    #######################################
    # Verify that no precision was lost during the computation
    assert verify.o(norm_xiL2_comp)

    verified = norm_xiL2_comp < norm_xiL2_save
    return verified, lemma_label, bounds


# In[ ]:


def norm_xi_fvap_sq(): # (Lemma 3.10 and 3.15)
    """
    Validate bounds associated with the approximate eigenfunction used in the
    stability argument.

    This function verifies three estimates involving the approximate Fourier
    eigenvector `fvap` and its associated approximate eigenvalue `laap`:

        1. A bound for the squared L2 norm of `fvap`.
        2. Bounds for the squared H1 norm of `fvap`.
        3. A bound for the squared L2 norm of the residual

               xiap = L_{the,twap}(fvap) - laap * fvap.

    The residual `xiap` measures how well `fvap` approximates an eigenvector
    of the linearized operator.  Its L2 norm provides the defect of the
    approximation.

    The computed quantities are enclosed are compared
    against the validated bounds stored in `BOUNDS`:

        ||fvap||_L2^2  < BOUNDS["fvap_stab_L2_sq"]
        ||fvap||_H1^2  < BOUNDS["fvap_stab_H1_sq"]
        ||fvap||_H1    > BOUNDS["rad_stab"]
        ||xiap||_L2^2  < BOUNDS["resi_stab_L2_sq"]

    """

    #######################################
    # Output information
    lemma_label = "Lemma 3.10 and 3.15"
    bound_labels = ['fvap_stab_L2_sq', 'fvap_stab_H1_sq', 'resi_stab_L2_sq', 'rad_stab']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]

    #######################################
    # Load data
    from explicit_funcs import fvap_der_and_Hil
    
    speed  = load_data.speed()
    coeffs = load_data.coeffs()
    laap   = load_data.eigen_value_approx()
    theta  = load_data.theta()
    fvap   = load_data.eigen_fv()
    
    fvap_L2sq_save = BOUNDS['fvap_stab_L2_sq']
    fvap_H1sq_save = BOUNDS['fvap_stab_H1_sq']
    resi_L2sq_save = BOUNDS['resi_stab_L2_sq']
    rad_stab       = BOUNDS["rad_stab"]
    
    #######################################
    # fvap_L2_sq
    fvap_L2sq_comp = auxiliar_funcs.norm_vector_sq(fvap)

    #######################################
    # fvap_H1_sq
    N = len(fvap) // 2
    assert N == len(coeffs)

    fvap_der, fvap_Hil = fvap_der_and_Hil(fvap, N)

    fvap_H1sq_comp = fvap_L2sq_comp + auxiliar_funcs.norm_vector_sq(fvap_der)
    
    #######################################
    # resi_L2_sq
    ## xiap = Lthelaap fvap - laap fvap,
    ## where:
    ## Lthelaap f = (c + twap) (f' + ithe f) + Hf - i f0
    ##                          ‾‾‾‾‾‾‾‾‾‾‾|_ this is what we call `aux_term`
    ##               ‾‾‾‾‾‾‾‾|_______________ this is what we call `ft`, we split it into two different parts:
    ##                                        the common indices with Hf and the rest 
    
    aux_term = fvap_der + I * theta * fvap
    
    # preallocate ft parts
    ft_common = speed * aux_term 
    ft_respos = vector(CBF, N)
    ft_resneg = vector(CBF, N)

    coeffs_half = [vj * ONE_DIV_2 for vj in coeffs]
    
    aux_term_p, aux_term_m = aux_term[:N], aux_term[N:]
    assert len(aux_term_p) == len(aux_term_m)
    
    for jj, vj in enumerate(coeffs_half, start=1):
        for ll, (coef_p, coef_m) in enumerate(zip(aux_term_p, aux_term_m)):            
            new_coef_p = vj * coef_p
            new_coef_m = vj * coef_m

            i_s = ll + jj
            i_d = ll - jj
            
            if i_s < N:
                ft_common[i_s    ] += new_coef_p
                ft_common[i_s + N] += new_coef_m
            else:
                ft_respos[i_s - N] += new_coef_p
                ft_resneg[i_s - N] += new_coef_m

            if i_d >= 0:
                ft_common[i_d    ] += new_coef_p
                ft_common[i_d + N] += new_coef_m
            else:
                i_reflect = - (i_d + 1)
                ft_common[i_reflect + N] += new_coef_p
                ft_common[i_reflect    ] += new_coef_m

    # `first_term` + Hf - i f0 common modes
    common_part_L = ft_common + fvap_Hil
    common_part_L[0] -= I * fvap[0]

    # Lf - lf 
    common_part = common_part_L - laap * fvap

    # norms
    norm_xi_sq_common = auxiliar_funcs.norm_vector_sq(common_part)
    norm_xi_sq_respos = auxiliar_funcs.norm_vector_sq(ft_respos)
    norm_xi_sq_resneg = auxiliar_funcs.norm_vector_sq(ft_resneg)
    
    norm_xisq_comp = norm_xi_sq_common + norm_xi_sq_respos + norm_xi_sq_resneg
    
    #######################################
    # Verifies that no precision was lost during the computation
    assert verify.o(fvap_L2sq_comp)
    assert verify.o(fvap_H1sq_comp)
    assert verify.o(norm_xisq_comp)
    
    cond_fvap_L2_sq = fvap_L2sq_comp < fvap_L2sq_save
    cond_fvap_H1_sq = fvap_H1sq_comp < fvap_H1sq_save
    cond_resi_L2_sq = norm_xisq_comp < resi_L2sq_save
    cond_lemma_3_15 = fvap_H1sq_comp > rad_stab**2
    verified = cond_fvap_L2_sq and cond_fvap_H1_sq and cond_resi_L2_sq and cond_lemma_3_15
    return verified, lemma_label, bounds


# In[ ]:





# In[ ]:


def coeffs_are_not_zero(): # (Lemma 4.3)
    """
    Checks that every cosine coefficient of the approximated traveling wave is never zero.
    """ 
    
    #######################################
    # Output information
    lemma_label = "Lemma 4.3"
    bounds = None

    #######################################
    # Load data
    coeffs = load_data.coeffs()
    
    #########################################
    # Check if every ball does not contain 0
    verified = all(not coef.contains_zero() for coef in coeffs)
    return verified, lemma_label, bounds


# In[ ]:


def enclose_roots(): # (Lemma C.3)
    """
    Verify validated enclosures for the roots of the polynomial associated with
    the approximate traveling wave.

    Let N = len(coeffs), and define the polynomial

        P(z) = z^N * (speed + sum_{j=1}^N (z^j + z^{-j}) * coeffs_j / 2).

    Using precomputed approximate roots `roots` and a common radius
    `radii_roots`, this function checks the following:

        1. For each approximate root `root` in `roots`, the closed disc
           D(root, radii_roots) is validated to contain a true root of P.

        2. Every such disc is contained in the open unit disc.

        3. All such discs are pairwise disjoint.

    To prove the root enclosure in (1), the function constructs coefficient-based
    bounds for the derivatives of P on the unit disc, namely numbers
    `big_bounds[d]` satisfying

        sup_{|z| <= 1} |dz^k P(z)| <= big_bounds[k].

    These bounds are then used by `root_polynomial_is_enclosed`.
    """

    #######################################
    # Output information
    lemma_label = "Lemma C.3"
    bound_labels = ['radii_roots']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]

    #######################################
    # Load data
    from explicit_funcs import palipoly_coeffs_constructor, big_bounds_compute_fast, root_polynomial_is_enclosed
    
    speed  = load_data.speed()
    coeffs = load_data.coeffs()
    roots  = load_data.polynomial_zeros()
    
    radii_roots = BOUNDS['radii_roots']

    #######################################
    # Construct the polynomial P associated to the traveling wave
    rk = len(coeffs) # N
    polynom = palipoly_coeffs_constructor(speed, coeffs)

    #######################################
    # Precompute auxiliary bounds for derivatives of P on |z|<=1.
    #
    # We build numbers `big_bounds[k]` such that, for each derivative order k,
    #     sup_{|z|<=1} |dz^k P(z)| <= big_bounds[k].
    #
    # For |z|<=1 and P(z)=sum_j a_j z^j, we have
    #     |dx^k P(z)| <= sum_{j>=k} |a_j| * j*(j-1)*...*(j-k+1),
    # hence the coefficient-based bound below.
    big_bounds = big_bounds_compute_fast(polynom)

    
    #######################################
    # Proof loop:
    # - condA: every proposed disc contains a root (validated via `root_polynomial_is_enclosed`)
    # - circ_dist: minimal distance margin to the unit circle (positive => all discs inside |z|<1)
    # - pair_dist: minimal separation margin between discs (positive => discs are disjoint)
    
    circ_dist = ONE
    pair_dist = TWO

    condA = True
    enclo_roots = [root.add_error(radii_roots) for root in roots]
    list_roots = list(zip(roots, enclo_roots))
    while list_roots and condA:
        # Take one approximate root and thicken it into a closed disc of radius `radii_roots`
        root, root_enclo = list_roots.pop()

        ###################################
        # (A) Enclosure: prove that P has a true root inside D(root, radii_roots)
        condA = root_polynomial_is_enclosed(root, polynom, radii_roots, big_bounds)
            
        ###################################
        # (B) Contained in the unit disc:
        # The margin to the unit circle is 1 - |D|, where |D| is the maximal modulus in the disc.
        len_dis = 1 - root_enclo.abs()
        circ_dist = circ_dist.min(len_dis)
        
        ###################################
        # (C) Disjoint discs: ensure D(root, r) does not intersect any remaining disc
        # We track the minimal margin of separation across all pairs.
        for _, other_root_enclo in list_roots:
            len_pair = (root_enclo-other_root_enclo).abs()
            pair_dist = pair_dist.min(len_pair)

        counter = rk - len(list_roots)
        print_iter_N(counter, text=f"{len(list_roots)}")
        

    #######################################
    # Final conditions:
    # - circ_dist > 0  => every disc is strictly inside the unit disc
    # - pair_dist > 0  => discs are pairwise disjoint
    condB = circ_dist > ZERO
    condC = pair_dist > ZERO

    assert verify.o(circ_dist)
    assert verify.o(pair_dist)
    
    verified = condA and condB and condC
    return verified, lemma_label, bounds
       


# In[ ]:


def proots_exis(): # (Lemma 4.12)
    """
    Verify a compatibility condition for the validated root enclosures of the
    polynomial associated with the approximate traveling wave.

    Let P be the polynomial built from `(speed, coeffs)`, let `roots` be the
    list of approximate roots, and let `radii_roots` be the common enclosure
    radius. For each approximate root, the function considers the corresponding
    enclosed disc and evaluates the interval quantity

        image_root(z) = - z^(N-1) / P'(z),
        
    where N is the number of approximate roots and z ranges over the enclosed
    root disc.

    For each enclosed root, the function checks that the enclosure of image_root cannot contain a real integer.
    More precisely, the condition is considered verified if at
    least one of the following holds:

        1. the imaginary part of `image_root` does not contain 0;
        2. the real part of `image_root` does not contain any integer.

    Therefore, the verification fails only if, for some enclosed root, the
    imaginary part may vanish and the real part may contain an integer.
    """

    #######################################
    # Output information
    from explicit_funcs import palipoly_coeffs_constructor, polynom_eval
    lemma_label = "Lemma 4.12"
    bound_labels = ['radii_roots']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]
    
    #######################################
    # Load data
    speed  = load_data.speed()
    coeffs = load_data.coeffs()
    roots  = load_data.polynomial_zeros()
    radii_roots = BOUNDS['radii_roots']
    
    #######################################
    # Construct the polynomial P associated to the traveling wave
    polynom = palipoly_coeffs_constructor(speed, coeffs)
   
    #######################################
    # Enclose each approximate root into a complex ball of radius `radii_roots`
    enclosed_roots = [root.add_error(radii_roots) for root in roots]

    #######################################
    # Lemma check over all enclosed roots
    #
    # For each enclosed root z, form:
    #     image_root = - z^(N-1) / P'(z).
    # We then test that it cannot be an integer real number under enclosure arithmetic.
    verified = True 
    while enclosed_roots and verified:
        root = enclosed_roots.pop()

        # Evaluate derivative P'(root) (idx_der=1) and build the enclosed image
        image_root = - root**(len(roots)-1) / polynom_eval(root, polynom, idx_der=1)

        # Separate real and imaginary parts (both are interval/ball quantities)
        image_root_real=image_root.real()
        image_root_imag=image_root.imag()

        # Pass condition:
        #  - either Im(image_root) is provably nonzero , OR
        #  - Re(image_root) is provably not an integer.
        verified = not (image_root_imag.contains_zero() and image_root_real.contains_integer())

    return verified, lemma_label, bounds


# In[ ]:


def proots_stab(): # (Lemma 5.25)
    """
    Verify a nondegeneracy condition for the validated root enclosures of the
    polynomial associated with the approximate traveling wave.

    Let P be the polynomial built from `(speed, coeffs)`, let `roots` be the
    list of approximate roots, let `radii_roots` be the common enclosure
    radius, and let `lamb` be the approximate eigenvalue.

    For each approximate root, the function considers the corresponding
    enclosed disc and evaluates the interval quantity

        aux_term(z) = z^(N-1) / P'(z),

    where N is the number of approximate roots and z ranges over the enclosed
    root disc.

    Using this quantity, the function forms the two images

        image_root1(z) = ( i*lamb - 1) * aux_term(z),
        image_root2(z) = (-i*lamb - 1) * aux_term(z).

    For each enclosed root, both images must avoid being compatible with a
    real integer. More precisely, an image passes the test if at least one of
    the following holds:

        1. its imaginary part does not contain 0;
        2. its real part does not contain any integer.

    Thus, failure occurs only if, for some enclosed root, the enclosure of an
    image may be real and its real part may contain an integer.
    """


    #######################################
    # Output information
    lemma_label = "Lemma 5.25"
    bound_labels = ['radii_roots']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]
    
    #######################################
    # Load data
    from explicit_funcs import palipoly_coeffs_constructor, polynom_eval
    
    speed  = load_data.speed()
    coeffs = load_data.coeffs()
    roots  = load_data.polynomial_zeros()
    lamb   = load_data.eigen_value_approx()
    
    radii_roots = BOUNDS['radii_roots']
    
    #######################################
    # Construct the polynomial P associated to the traveling wave
    polynom = palipoly_coeffs_constructor(speed, coeffs)

    #######################################
    # Enclose each approximate root into a complex ball of radius `radii_roots`
    enclosed_roots = [root.add_error(radii_roots) for root in roots]

    #######################################
    # Lemma check over all enclosed roots
    #
    # For each enclosed root z we compute aux_term(z) = z^(N-1)/P'(z) and test the two
    # scaled images (λ-1)*aux_term and (-λ-1)*aux_term against the "integer real" obstruction.
    verified = True 
    while enclosed_roots and verified:
        root = enclosed_roots.pop()
        
        # Compute aux_term(z) = z^(N-1) / P'(z) with enclosure arithmetic
        aux_term = root**(len(roots)-1) / polynom_eval(root, polynom, idx_der=1)

        # Two images: (±λ - 1) * aux_term, for each sign of lamb.
        image_root1 = ( I*lamb - 1) * aux_term
        image_root2 = (-I*lamb - 1) * aux_term

        # Separate real/imaginary parts (both are interval/ball quantities)
        image_root_real1=image_root1.real()
        image_root_imag1=image_root1.imag()
        image_root_real2=image_root2.real()
        image_root_imag2=image_root2.imag()

        # Pass condition for each image:
        #  - either Im(image_root) is provably nonzero (does NOT contain 0 exactly), OR
        #  - Re(image_root) is provably not an integer (does NOT contain any integer).
        #
        # Failure happens only if it could be exactly real AND could hit an integer.
        veri1 = not (image_root_imag1.contains_zero() and image_root_real1.contains_integer())
        veri2 = not (image_root_imag2.contains_zero() and image_root_real2.contains_integer())

        # We require both images to satisfy the condition for every root
        verified = veri1 and veri2

    return verified, lemma_label, bounds


# In[ ]:





# In[ ]:


def beta_max(): # (Lemma 4.18)
    """
    Validate the uniform bound for
    
        beta(x) = 1 / (c + v(x))
    
    on [0,pi] and beta(x) > 0.
    """
    order = 2

    #######################################
    # Output information
    lemma_label = "Lemma 4.18"
    bound_labels = ['beta_max']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]
    
    #######################################
    # Load data
    from explicit_funcs import beta_func_constructor
    
    speed  = load_data.speed()
    coeffs = load_data.coeffs()
    
    beta_max_save = BOUNDS['beta_max']

    #######################################
    # Construct beta(x) = 1 / (c + v(x)) together with its derivatives up to the prescribed order.
    beta_func = beta_func_constructor(speed, coeffs, order=order)

    #######################################
    # Domain construction
    domain = [ZERO, PI]

    # Choose a kth-order enclosure method for the supremum computation
    method = image_1D.taylor_method(order)

    #######################################
    # Verify the supremum bound using an adaptive enclosure algorithm
    satisf_boun = auxiliar_funcs.max_prev_esti(domain,         beta_func, beta_max_save, method)
    is_positive = auxiliar_funcs.max_prev_esti(domain, RBF(-1)*beta_func,          ZERO, method)
    
    verified = satisf_boun and is_positive
    return verified, lemma_label, bounds


# In[ ]:


def kappa1(): # (Lemma 4.19)
    """
    Verifies an L1 bound for
        kappa1(x) := int_0^x beta(y)^2 dy,
    where
        beta(x) := 1 / (c + v(x)),
    and v is an approximation of the traveling wave with speed c.
    """
    abs_tol = RBF('1e-4') 
    rel_tol = RBF('1e-5')

    #######################################
    # Output information
    lemma_label = "Lemma 4.19"
    bound_labels = ['kappa1']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]

    #######################################
    # Load data
    from explicit_funcs import aux_kappa1_func_constructor
    
    speed  = load_data.speed()
    coeffs = load_data.coeffs()
    
    kappa1_save = BOUNDS['kappa1']

    #######################################
    # Construct betamod(x) = (pi - x) * beta(x)^2 (and derivatives if needed by the method)
    order = 4
    aux_kappa1_func = aux_kappa1_func_constructor(speed, coeffs, order)
    
    #######################################
    # Domain construction
    domain = [ZERO, PI]

    # Choose a kth-order enclosure method for the supremum computation
    method = integral_1D.gauss_lege_2dots()
    
    #######################################
    kappa1_comp = auxiliar_funcs.integ_adaptive_1D(domain, aux_kappa1_func, method, abs_tol, rel_tol)

    #######################################
    # Verifies that no precision was lost during the computation
    assert verify.o(kappa1_comp)
    
    verified = kappa1_comp < kappa1_save
    return verified, lemma_label, bounds


# In[ ]:


def betamod_L2sq(): # (Lemma 4.20)
    """
    Validate the L^2 bound for
    
    beta'(x) + i beta(x)^2.
    """
    abs_tol = RBF('1e-2') 
    rel_tol = RBF('1e-4')

    #######################################
    # Output data
    lemma_label = "Lemma 4.20"
    bound_labels = ['beta_mod']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]
    
    #######################################
    # Load data
    from explicit_funcs import betamod_func_sq_constructor
    
    speed  = load_data.speed()
    coeffs = load_data.coeffs()
    
    betamod_L2sq_save = BOUNDS['beta_mod']

    #######################################
    # Construct integrand
    order = 4
    integrand = betamod_func_sq_constructor(speed, coeffs, order)

    # Construct domain
    domain = [ZERO, PI]

    # Enclosure method for the integral (adaptive Gauss–Legendre)
    method = integral_1D.gauss_lege_2dots()

    # Compute a rigorous enclosure of the integral
    betamod_L2sq_comp = 2 * auxiliar_funcs.integ_adaptive_1D(domain, integrand, method, abs_tol, rel_tol)

    #######################################
    # Verifies that no precision was lost during the computation
    assert verify.o(betamod_L2sq_comp)
    
    verified = betamod_L2sq_comp < betamod_L2sq_save
    return verified, lemma_label, bounds


# In[ ]:


def kappa2(): # (Lemma 4.21)
    """
    Validate the L1 bound for kappa2.
    
    The quantity is computed through the triangular integral obtained after
    exchanging the order of integration.
    """
    abs_tol = RBF("1e-1")

    #######################################
    # Output information
    lemma_label = "Lemma 4.21"
    bound_labels = ['kappa2']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]

    #######################################
    # Load data
    from explicit_funcs import beta_sq_func_constructor, betamod_func_sq_constructor
    
    speed  = load_data.speed()
    coeffs = load_data.coeffs()
    
    kappa2_save = BOUNDS['kappa2']

    #######################################
    order = 4
    beta_sq = beta_sq_func_constructor(speed, coeffs, order)

    # Nonnegative diota := |beta'(x) + i*beta(x)^2|^2 = beta'(x)^2 + beta(x)^4
    diota = betamod_func_sq_constructor(speed, coeffs, order)

    # iota(x) = int_0^x diota(y) dy
    # betaiota_L1 = int_0^pi beta(y)^2 (iota(pi) - iota(y)) dy
    kappa2_comp = auxiliar_funcs.triangular_integral(diota, beta_sq, abs_tol)

    # Non precision lost
    assert verify.o(kappa2_comp)

    verified = kappa2_comp < kappa2_save
    return verified, lemma_label, bounds


# In[ ]:


def kappa3(): # (Lemma 5.27)
    """
    Validate the L1 bound for kappa3.
    """
    
    abs_tol = RBF('1e-4') 
    rel_tol = RBF('1e-4')

    #######################################
    # Output data
    lemma_label = "Lemma 5.27"
    bound_labels = ['kappa3','kappa1']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]

    #######################################
    # Load data
    from explicit_funcs import beta_func_constructor
    
    speed  = load_data.speed()
    coeffs = load_data.coeffs()
    laap   = load_data.eigen_value_approx()
    
    kappa1 = BOUNDS['kappa1']
    kappa3_save  = BOUNDS['kappa3']
    
    #######################################
    # Enclose int_0^pi beta(y) y
    order = 4
    beta_func = beta_func_constructor(speed, coeffs, order)

    domain = [ZERO, PI]

    method = integral_1D.gauss_lege_2dots()

    # Compute a rigorous enclosure of the integral
    int_0pi_be = auxiliar_funcs.integ_adaptive_1D(domain, beta_func, method, abs_tol, rel_tol)

    #######################################
    # Compute the upper multiplicative bound
    ## C = e^(2 |Re(laap)| int_0^pi beta(y)^2 dy) 
    C = (TWO * laap.real().abs() * int_0pi_be).exp()

    #######################################
    # Construct result and verify the statement
    kappa3_comp = C * kappa1
    assert verify.o(kappa3_comp)

    verified = kappa3_comp < kappa3_save
    return verified, lemma_label, bounds


# In[ ]:


def kappa4(): # (Lemma 5.28)
    """
    Validate the L1 bound for the kappa4 functions.
    """
    abs_tol_1 = RBF('1e-4') 
    rel_tol_1 = RBF('1e-4')
    abs_tol_2 = RBF("1e-2")

    #######################################
    # Output data
    lemma_label = "Lemma 5.28"
    bound_labels = ['kappa4p','kappa4m','kappa1']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]

    #######################################
    # Load data
    from explicit_funcs import beta_func_constructor, betamod_laap_sq_func_constructor, beta_sq_func_constructor
    
    speed  = load_data.speed()
    coeffs = load_data.coeffs()
    laap   = load_data.eigen_value_approx()
    theta  = load_data.theta()
    
    kappa1 = BOUNDS['kappa1']
    kappa4p_save = BOUNDS['kappa4p']
    kappa4m_save = BOUNDS['kappa4m']


    #######################################
    # Enclose int_0^pi beta(y) y
    order = 4
    beta_func = beta_func_constructor(speed, coeffs, order)

    domain = [ZERO, PI]
    method = integral_1D.gauss_lege_2dots()

    # Compute a rigorous enclosure of the integral
    int_0pi_besq = auxiliar_funcs.integ_adaptive_1D(domain, beta_func, method, abs_tol_1, rel_tol_1)

    #######################################
    # Compute the upper multiplicative bound
    ## C = e^(2 |Re(laap)| int_0^pi beta(y)^2 dy) 
    C = (TWO * laap.real().abs() * int_0pi_besq).exp()

    #######################################
    #
    integ_p, integ_m = betamod_laap_sq_func_constructor(speed, coeffs, laap, theta, order)
    beta_sq = beta_sq_func_constructor(speed, coeffs, order)

    kappa4p_comp = C * auxiliar_funcs.triangular_integral(integ_p, beta_sq, abs_tol_2)
    kappa4m_comp = C * auxiliar_funcs.triangular_integral(integ_m, beta_sq, abs_tol_2)
    
    # Non precision lost
    assert verify.o(kappa4p_comp)
    assert verify.o(kappa4m_comp)

    verified = (kappa4p_comp < kappa4p_save) and (kappa4m_comp < kappa4m_save)
    return verified, lemma_label, bounds
    


# In[ ]:





# In[ ]:


def hk_norm_exis_H2sq(): # (Lemma 4.25)
    """
    Validate the bound for the summed L^2 second derivatives of the hk functions (4.44). 
    """
    #######################################
    # Output information
    lemma_label = "Lemma 4.25"
    bound_labels = ['beta_max','beta_mod','kappa2','hkH2']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]
    
    #######################################
    # Load data
    from explicit_funcs import compute_hk_norms_exis_ap, beta_dx_func_sq_constructor, real_gk_funcs_constructor

    speed  = load_data.speed()
    coeffs = load_data.coeffs()
    gkLnfr = load_data.Linf_gk_bounds()
    
    betamod_L2sq = BOUNDS['beta_mod']
    betamax      = BOUNDS['beta_max']
    kappa2       = BOUNDS['kappa2']
    hkH2sq_save  = BOUNDS['hkH2']

    #######################################

    domain = [ZERO, PI]
    order = 2
    method = image_1D.taylor_method(order)
    half = ONE_DIV_2

    real_gk_funcs = real_gk_funcs_constructor(coeffs, order)
    assert len(real_gk_funcs) == len(gkLnfr)
    for kk, (func, bound) in enumerate(zip(real_gk_funcs, gkLnfr), start=1):
        cond1 = auxiliar_funcs.max_prev_esti(domain,        func, bound, method, verbose=0)
        cond2 = auxiliar_funcs.max_prev_esti(domain, (-ONE)*func, bound, method, verbose=0)

        print_iter_N(kk, text=f'{kk}')

        if not cond1 or not cond2:
            raise ValueError(
                f"Function g{kk} fails: {cond1} and {cond2}"
            )
        
    #######################################

    N=len(coeffs)
    acc_sq_norm = betamod_L2sq
    I1ac = ZERO; I2ac = ZERO; I3ac = ZERO;
    for kk, vk in enumerate(coeffs,start=1):
        # --- I1 ---
        dgk_L2sq = ZERO
        for mm, vkm in enumerate(coeffs[kk:], start=1):
            dgk_L2sq += (mm*vkm)**2
        I1sq = PI * dgk_L2sq * betamax**2 * kk**2

        # --- I2 ---
        re_part = gkLnfr[kk-1]
        im_part = kk * vk.abs()*half
        sq_part = re_part**2 + im_part**2
        I2sq = sq_part * betamod_L2sq

        # --- I3 ---
        gkk_L2sq = half * vk**2
        for mm, vkm in enumerate(coeffs[kk:], start=1):
            gkk_L2sq  += vkm**2
        I3sq = PI * kappa2 * gkk_L2sq * kk**2
        
        # ----------
        II = [I1sq, I2sq, I3sq]
        I1ac += I1sq; I2ac += I2sq; I3ac += I3sq;
        acc_sq_norm += sum(x**half for x in II)**2

    assert verify.o(acc_sq_norm)
    
    verified = acc_sq_norm < hkH2sq_save
    return verified, lemma_label, bounds


# In[11]:


def hk_norm_stab_H1sq(): # (Lemma 5.31)
    """
    Validate the L^2 and `H^1` bounds for the singular stability modes h_k^+ and h_k^- (5.28).
    """
    #######################################
    # Output information
    lemma_label = "Lemma 5.31"
    bound_labels = ['kappa3','kappa4p','kappa4m','hkL2p','hkL2m','hkH1p','hkH1m']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]
    
    #######################################
    # Load data
    from explicit_funcs import compute_hk_norms_stab_ap
    
    ode_values  = load_data.ode_stab_sing_values()
    ode_bounds0 = load_data.ode_stab_sing_bounds()
    
    JL2L2sq  = BOUNDS['kappa3']
    JL2H1psq = BOUNDS['J1p']
    JL2H1msq = BOUNDS['J1m']
    hkL2psq_save = BOUNDS['hkL2p']
    hkL2msq_save = BOUNDS['hkL2m']
    hkH1psq_save = BOUNDS['hkH1p']
    hkH1msq_save = BOUNDS['hkH1m']

   
    def compute_hk_norms(hknorms_sq, residue_sq):
        """ We add the residual part """
        half = ONE_DIV_2
        assert len(hknorms_sq) == len(residue_sq)
        return sum( (a**half + b**half)**2 for a, b in zip(hknorms_sq, residue_sq))

    def aux_resi_sort(residues):
        """ Sort the residues """
        n = len(residues) // 2
        nm1 = n - 1

        residues_p = vector(RBF, n)
        residues_m = vector(RBF, n)
    
        residues_p[:nm1] = residues[:nm1]
        residues_m[:nm1] = residues[nm1:2*nm1]
        residues_p[-1] = residues[-2]
        residues_m[-1] = residues[-1]

        return residues_p, residues_m

    #######################################
    hk_norm_ap_p0, hk_norm_ap_m0, hk_norm_ap_p1, hk_norm_ap_m1 = compute_hk_norms_stab_ap(ode_values)
    
    residues_p0, residues_m0 = aux_resi_sort(ode_bounds0)
    residues_p1 = residues_p0 * JL2H1psq / JL2L2sq
    residues_m1 = residues_m0 * JL2H1msq / JL2L2sq

    hkL2psq_comp = compute_hk_norms(hk_norm_ap_p0, residues_p0)
    hkL2msq_comp = compute_hk_norms(hk_norm_ap_m0, residues_m0)
    hkH1psq_comp = compute_hk_norms(hk_norm_ap_p1, residues_p1)
    hkH1msq_comp = compute_hk_norms(hk_norm_ap_m1, residues_m1)
    
    # Non precision lost
    assert verify.o(hkL2psq_comp)
    assert verify.o(hkL2msq_comp)
    assert verify.o(hkH1psq_comp)
    assert verify.o(hkH1msq_comp)

    ver0 = (hkL2psq_comp < hkL2psq_save) and (hkL2msq_comp < hkL2msq_save)
    ver1 = (hkH1psq_comp < hkH1psq_save) and (hkH1msq_comp < hkH1msq_save)
    verified = ver0 and ver1
    return verified, lemma_label, bounds
    


# In[ ]:





# In[ ]:


def svd_exis(): # (Lemma 4.16)
    """
    Validate the lower bound for the smallest singular value of M_exis^tor.
    """
    #######################################
    # Output data
    lemma_label = "Lemma 4.16"
    bound_labels = ['svd1_exis']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]
    
    #######################################
    # Load data
    from explicit_funcs import construct_exis_mat, mat_exis_add_radii
    
    ode_values = load_data.ode_exis_values()
    ode_bounds = load_data.ode_exis_bounds()
    step_mat_V = load_data.ode_exis_stepmat()
    min_svd_save = BOUNDS['svd1_exis']

    #######################################
    # Check structural invariants
    n = len(ode_bounds)
    assert verify.data_block(ode_values, n)
    
    #######################################
    # Construct the thin existence matrix
    ord_mat = n
    mat_thin = construct_exis_mat(ode_values)
    
    #######################################
    # Dimension checks (avoid opaque Sage errors later)
    assert verify.square_mat_order(mat_thin, ord_mat, 'mat_exis')
    assert verify.square_mat_order(step_mat_V, ord_mat, 'V_exis')

    #######################################
    # Compute perturbation bound due to the approximated functions `hk`
    perturb = mat_exis_add_radii(mat_thin, ode_bounds)
    
    #######################################
    # Build V* M* M V
    aux_diag = mat_thin * step_mat_V
    ti_S = aux_diag.conjugate().transpose() * aux_diag
    M = mat_thin.conjugate().transpose() * mat_thin
    S = auxiliar_funcs.S_from_ti_S_for_Gersh(ti_S, M, step_mat_V)

    #######################################
    # Gershgorin lower bounds for a Hermitian matrix:
    #   λ_min >= min_j ( a_jj - sum_{k!=j} |a_jk| )
    min_svd_comp = auxiliar_funcs.gersh_smallest_eig(S)

    #######################################
    # Subtract perturbation bound and verify strict inequality
    min_svd_comp -= perturb
    
    assert verify.o(min_svd_comp)
    
    verified = min_svd_save < min_svd_comp
    return verified, lemma_label, bounds


# In[ ]:


def svd_regu(): # (Lemma 5.16)
    """
    Validate the positivity for the smallest singular value of M_stab^tor(1.035i).
    """
    #######################################
    # Output data
    lemma_label = "Lemma 5.16"
    bounds = None
    
    #######################################
    # Load data
    from explicit_funcs import construct_stab_mat, mat_stab_add_radii
    
    ode_values = load_data.ode_stab_regu_values() # must be: pp,pm,mp,mm,ph,mh
    ode_bounds = load_data.ode_stab_regu_bounds()
    step_mat_V = load_data.ode_stab_regu_stepmat()

    #######################################
    # Check structural invariants
    n = len(ode_bounds)
    assert verify.data_block(ode_values, n)

    #######################################
    # Construct the thin stability matrix
    mat_thin = construct_stab_mat(ode_values)
    ord_mat = 1 + n//2

    #######################################
    # Dimension checks (avoid opaque Sage errors later)
    assert verify.square_mat_order(mat_thin, ord_mat, 'mat_stab')
    assert verify.square_mat_order(step_mat_V, ord_mat, 'V_exis')

    #######################################
    # quantity from the ODE residues
    perturb = mat_stab_add_radii(mat_thin, ode_bounds)            

    #######################################
    # Build V* M* M V
    aux_diag = mat_thin * step_mat_V
    ti_S = aux_diag.conjugate().transpose()*aux_diag 
    M    = mat_thin.conjugate().transpose() * mat_thin
    S    = auxiliar_funcs.S_from_ti_S_for_Gersh(ti_S, M, step_mat_V)
    
    min_val = auxiliar_funcs.gersh_smallest_eig(S)

    #######################################
    # Subtract perturbation bound and verify strict inequality
    min_val -= perturb
    
    assert verify.o(min_val)
    
    verified = min_val > ZERO
    return verified, lemma_label, bounds
    


# In[ ]:


def svd_sing(): # (Lemma 5.22)
    """
    Let M = M_stab^tor(laap)
    It proves an upper bound for the first (smallest) singular value of M and a lower bound for
    the second one.
    """
    #######################################
    # Output data
    lemma_label = "Lemma 5.22"
    bound_labels = ['svd_stab1','svd_stab2']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]
    
    #######################################
    # Load data
    from explicit_funcs import construct_stab_mat, mat_stab_add_radii
    
    ode_values = load_data.ode_stab_sing_values()
    ode_bounds = load_data.ode_stab_sing_bounds()
    step_mat_V = load_data.ode_stab_sing_stepmat()
    
    svd1_save = BOUNDS['svd_stab1']
    svd2_save = BOUNDS['svd_stab2']

    #######################################
    #  Check structural invariants
    n = len(ode_bounds)
    assert verify.data_block(ode_values, n)

    #######################################
    # Construct the thin stability matrix
    mat_thin = construct_stab_mat(ode_values)
    ord_mat = 1 + n//2

    #######################################
    # Dimension checks (avoid opaque Sage errors later)
    assert verify.square_mat_order(mat_thin, ord_mat, 'mat_stab_sing')
    assert verify.square_mat_order(step_mat_V, ord_mat, 'V_exis')

    #######################################
    # quantity from the ODE residues
    perturb = mat_stab_add_radii(mat_thin, ode_bounds)            

    #######################################
    # Build V* M* M V
    aux_diag = mat_thin * step_mat_V
    ti_S = aux_diag.conjugate().transpose()*aux_diag 
    M = mat_thin.conjugate().transpose() * mat_thin
    S = auxiliar_funcs.S_from_ti_S_for_Gersh(ti_S, M, step_mat_V)
    
    svd1_comp, svd2_comp = auxiliar_funcs.gersh_gap_one_two(S)

    svd1_comp += perturb
    svd2_comp -= perturb
    assert verify.o(svd1_comp)
    assert verify.o(svd2_comp)
    
    verified = (svd1_comp < svd1_save) and (svd2_comp > svd2_save)
    return verified, lemma_label, bounds
    


# In[ ]:


def prod_fv_u1(): # (Lemma 5.23)
    """
    Verify a lower bound estimate for a dot product.
    """

    #######################################
    # Output information
    lemma_label = "Lemma 5.23"
    bound_labels = ['svd_stab1', 'svd_stab2', 'Jfvu1', 'kappa3', 'fvap_stab_L2_sq']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]
    
    #######################################
    # Load data
    from explicit_funcs import construct_stab_mat, bfv_constructor, div_twopi_sqrt_vec, col_norms, aux_sort_stab_mat
    ode_values = load_data.ode_stab_sing_values()
    ode_bounds = load_data.ode_stab_sing_bounds()
    Jfv_values = load_data.ode_stab_Jfv_values()
    Jfv_bounds = load_data.ode_stab_Jfv_bounds()
    u1_ap      = load_data.eigen_u1() # with norm 1
    
    sv1_sq_upp = BOUNDS['svd_stab1']
    sv2_sq_low = BOUNDS['svd_stab2']
    kappa3     = BOUNDS['kappa3']
    fvap_L2sq  = BOUNDS['fvap_stab_L2_sq']
    fv_u1_save = BOUNDS['Jfvu1']

    #######################################
    # Construct the thin stability matrix
    mat_thin = construct_stab_mat(ode_values)
        
    mat_order = mat_thin.ncols()
    nm1 = mat_order//2 - 1
    nn = len(ode_bounds)
    assert nn // 4 == nm1
    ord_mat = 2*nm1 + 2

    #######################################
    # xi_u1 = xi_thin + xi_resi
    Mast_u1 = mat_thin.conjugate().transpose() * u1_ap
    a_sq = auxiliar_funcs.norm_vector_sq(Mast_u1)
    aux_vec_xi_thin = mat_thin * Mast_u1 - a_sq * u1_ap
    xi_thin = aux_vec_xi_thin.norm()

    regular_resii = aux_sort_stab_mat(ode_bounds)
    RR = div_twopi_sqrt_vec(regular_resii)
    MM = col_norms(mat_thin) 

    Mast_u1_abs = vector(RBF, [a.abs() for a in Mast_u1])
    
    MRu = auxiliar_funcs.prod_vector(MM, RR)
    RMu = auxiliar_funcs.prod_vector(Mast_u1_abs, RR)
    RRu = auxiliar_funcs.prod_vector(RR, RR)
    xi_resi = MRu + RMu + RRu
    
    xi_u1 = 2*(xi_thin + xi_resi)

    #######################################
    # Estimate delta
    gap = sv2_sq_low - sv1_sq_upp
    de = (2**ONE_DIV_2) * xi_u1 / gap

    #######################################
    # Construct bfv and compute the approximate product
    bfv = bfv_constructor(Jfv_values, nm1, ord_mat)
    ap_prod = auxiliar_funcs.prod_vector(bfv, u1_ap)
    
    #######################################
    # Compute bfvap := b(Fap+,Fap-)

    Jfv_resi_sq = sum(Jfv_bounds) / TWOPI 
    kappa3_fvap = (kappa3 * fvap_L2sq)**ONE_DIV_2
    
    fv_u1_comp = ap_prod.abs() - Jfv_resi_sq ** ONE_DIV_2 - kappa3_fvap * de

    assert verify.o(fv_u1_comp)

    verified = (fv_u1_comp > fv_u1_save) and (xi_u1 < gap)
    return verified, lemma_label, bounds
    


# In[ ]:





# In[ ]:


def ode_resis_exis(): # (Lemma D.1)
    """
    Verifies an L2 bound for the difference between the hk functions (4.44) and our explicit approximations.
    """

    #######################################
    # Output data
    lemma_label = "Lemma D.1"
    bound_labels = ['kappa1']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]
    
    #######################################
    # Load data
    from explicit_funcs import compute_resis_exis
    
    speed      = load_data.speed()
    coeffs     = load_data.coeffs()
    ode_values = load_data.ode_exis_values()
    ode_bounds = load_data.ode_exis_bounds()
    
    JL2L2sq = BOUNDS['kappa1']

    #######################################
    # Check structural invariants
    n = len(ode_bounds) 
    assert len(coeffs) + 1 == n
    assert verify.data_block(ode_values, n)

    #######################################
    # Compute residues
    residues_sq = compute_resis_exis(speed, coeffs, ode_values)

    # Compare against bounds with a strict inequality
    residues_sq *= JL2L2sq
    
    # Avoid zip-truncation and no prec lost
    assert verify.oo(residues_sq, len(ode_bounds))

    verified = all(uk < vk for uk, vk in zip(residues_sq, ode_bounds))
    return verified, lemma_label, bounds
    


# In[ ]:


def ode_resis_regu(packed=None): # (Lemma D.7)
    """
    Verifies an L2 bound for the difference between the hk functions (D.18) and our explicit approximations,
    with lambda = 1.035i.
    """

    #######################################
    # Output data
    lemma_label = "Lemma D.7"
    bound_labels = ['kappa1']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]

    #######################################
    # Load data
    from explicit_funcs import compute_resis_stab
    
    speed  = load_data.speed()
    coeffs = load_data.coeffs()
    theta  = load_data.theta()
    if packed is None:
        
        lamb       = load_data.eigen_value_regular()
        ode_values = load_data.ode_stab_regu_values() # must be: pp,pm,mp,mm,ph,mh
        ode_bounds = load_data.ode_stab_regu_bounds()
        
        JL2L2sq = BOUNDS['kappa1']
    else:
        lemma_label, bounds, lamb, ode_values, ode_bounds, JL2L2sq = packed

    #######################################
    # Check structural invariants
    n = len(ode_bounds)
    assert 4*len(coeffs) + 2 == n
    assert verify.data_block(ode_values, n)
    
    #######################################
    # Compute residues
    residues_sq = compute_resis_stab(lamb, theta, speed, coeffs, ode_values)
 
    # Compare against bounds with a strict inequality
    residues_sq *= JL2L2sq

    assert verify.oo(residues_sq, len(ode_bounds))
    
    verified = all(uk < vk for uk, vk in zip(residues_sq, ode_bounds))
    return verified, lemma_label, bounds
    




# In[ ]:


def ode_resis_sing(): # (Lemma D.8)
    """
    Verifies an L2 bound for the difference between the hk functions (D.18) and our explicit approximations,
    with lambda = laap.
    """

    lemma_label = "Lemma D.8"
    bound_labels = ['kappa3']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]
    
    lamb = load_data.eigen_value_approx()
    ode_values = load_data.ode_stab_sing_values()
    ode_bounds = load_data.ode_stab_sing_bounds()
    JL2L2sq = BOUNDS['kappa3']

    packed = [lemma_label, bounds, lamb, ode_values, ode_bounds, JL2L2sq]
    
    return ode_resis_regu(packed=packed)


# In[ ]:


def ode_resis_Jfv(): # (Lemma D.10)
    """
    See Section D.2.1.
    """

    #######################################
    # Output data
    lemma_label = "Lemma D.10"
    bound_labels = ['kappa3']
    bounds = [BOUNDS_STR[lab] for lab in bound_labels]

    #######################################
    # Load data
    from explicit_funcs import compute_resis_Jfv
    
    speed      = load_data.speed()
    coeffs     = load_data.coeffs()
    lamb       = load_data.eigen_value_approx()
    ode_values = load_data.ode_stab_Jfv_values()
    ode_bounds = load_data.ode_stab_Jfv_bounds()
    fvap       = load_data.eigen_fv()
    theta      = load_data.theta()
    
    JL2L2sq = BOUNDS['kappa3']

    #######################################
    # Check structural invariants
    n = len(ode_bounds)
    assert n == 2
    assert verify.data_block(ode_values, n + 2)

    #######################################
    # Compute the bounds
    resip, resim = compute_resis_Jfv(lamb, theta, fvap, speed, coeffs, ode_values)

    resip *= JL2L2sq
    resim *= JL2L2sq

    assert verify.o(resip)
    assert verify.o(resim)
    
    verified = (resip < ode_bounds[0]) and (resim < ode_bounds[1])
    return verified, lemma_label, bounds
    




# In[ ]:


residual_norms = [
    norm_xi_twap_sq, 
    norm_xi_fvap_sq,
]

non_degen_coeffs = [
    coeffs_are_not_zero, 
    enclose_roots, 
    proots_exis, 
    proots_stab,
]

integral_bounds = [
    beta_max,     
    kappa1,       
    betamod_L2sq, 
    kappa2, 
    kappa3,       
    kappa4,       
]

hk_norms = [
    hk_norm_exis_H2sq, 
    hk_norm_stab_H1sq,
]

svd_results = [
    svd_exis, 
    svd_regu, 
    svd_sing, 
    prod_fv_u1,
]

ode_resis_results = [
    ode_resis_exis, 
    ode_resis_regu, 
    ode_resis_sing, 
    ode_resis_Jfv,
]


# In[ ]:


def substituting_estimates():
    """
    Checks that every substitution in the paper is numerically consistent.
    Returns (verified, description) where description is None if all checks pass,
    or a string listing the failed checks otherwise.
    """
    checks = []   # list of (bool, str) pairs: (passed, description)

    def add_check(condition, description):
        checks.append((condition, description))

    zet2 = PI**2 / RBF(6)
    theta = load_data.theta()
    half = ONE_DIV_2

    #######################################
    # --- Section 2, Proposition 2.5 ---
    rad_exis = BOUNDS['rad_exis']
    ell_exis = BOUNDS['lip_exis']
    
    qua = (2*zet2)**half
    lin = BOUNDS['Lexis_inv']
    aff = BOUNDS['resi_exis_L2']
    
    qq = qua * lin
    aa = aff * lin

    discriminant = 1 - 4 * aa * qq

    minimal_rad = (1 - discriminant**half) / (2*qq)
    maximal_rad = (1 + discriminant**half) / (2*qq)

    xx = rad_exis

    cont_exis = aa + qq * xx**2
    cont_cond = cont_exis < xx
    
    lips_exis = 2 * qq * xx
    lips_cond = lips_exis < ell_exis

    verify.o(minimal_rad)
    verify.o(cont_exis)
    verify.o(lips_exis)

    add_check(discriminant > 0, f"Texis discriminant not positive: {discriminant}")
    add_check(cont_cond, f"Radius for Texis not valid: {cont_exis} not less than {xx}")
    add_check(lips_cond, f"Texis not Lipschitz: {lips_exis}")

    #######################################
    # --- Corollary 3.9 ---
    Dthe_save = BOUNDS['Dthe_op']
    
    coef = 2 * zet2**half * (1 + theta**2)**half
    Dthe_comp = coef * rad_exis
    assert verify.o(Dthe_comp)

    add_check(Dthe_comp < Dthe_save, f"Dthe norm not valid: {Dthe_comp} not less than {Dthe_save}")    

    #######################################
    # --- Lemma 4.23 ---
    C1 = BOUNDS['beta_max']
    C2 = (BOUNDS['beta_mod'] * PI / RBF(12))**half
    C3 = BOUNDS['kappa2']**half
    JH1H2_save = BOUNDS['JH1H2']

    JH1H2_comp = C1 + C2 + C3

    assert verify.o(JH1H2_comp)
    add_check(JH1H2_comp < JH1H2_save, f"Jexis H1-H2 norm not valid: {JH1H2_comp} not less than {JH1H2_save}")    
    
    #######################################
    # --- Section 4.2 ---

    kappa1_sq    = BOUNDS['kappa1']
    hkH2_sq      = BOUNDS['hkH2']
    svd1_exis_sq = BOUNDS['svd1_exis']
    JH1H2        = BOUNDS['JH1H2']
    L_inv_save   = BOUNDS['Lexis_inv']

    # C1^2 := JL2L2^2 / s1^2  sum norm(dx^2 hk)^2
    C1 = (kappa1_sq * hkH2_sq / svd1_exis_sq)**half
    
    # C2 := JH1H2
    C2 =  JH1H2

    # C3 := C1 / sqrt(2 pi) + C2
    aux_PI = (2 * PI)**half
    C3 = C1 / aux_PI + C2

    assert verify.o(C3)
    add_check(C3 < L_inv_save, f"L_inv_exis norm not valid: {C3} not less than {L_inv_save}")    

    #######################################
    # --- Lemma 5.30 ---
    kappa4p  = BOUNDS['kappa4p']
    kappa4m  = BOUNDS['kappa4m']
    beta_max = BOUNDS['beta_max']
    J1p_save = BOUNDS['J1p']
    J1m_save = BOUNDS['J1m']
    
    J1p_comp = kappa4p**half + beta_max
    J1m_comp = kappa4m**half + beta_max

    assert verify.o(J1p_comp)
    assert verify.o(J1m_comp)
    add_check(J1p_comp < J1p_save, f"Jp L2-H1 norm not valid: {J1p_comp} not less than {J1p_save}")  
    add_check(J1m_comp < J1m_save, f"Jm L2-H1 norm not valid: {J1m_comp} not less than {J1m_save}")  
    
    #######################################
    # --- Proposition 5.32 ---
    C0p = BOUNDS['hkL2p']
    C0m = BOUNDS['hkL2m']
    C1p = BOUNDS['hkH1p']
    C1m = BOUNDS['hkH1m']

    CJ0 = BOUNDS['kappa3']
    CJp = BOUNDS['J1p']
    CJm = BOUNDS['J1m']

    sv2 = BOUNDS['svd_stab2']

    aux_coef = (CJ0 / (sv2 * TWOPI))**half
    ti_C0p_00 = aux_coef * C0p**half
    ti_C0m_00 = aux_coef * C0m**half
    ti_C1p_00 = aux_coef * C1p**half
    ti_C1m_00 = aux_coef * C1m**half

    ti_C0p_0p = CJ0**half
    ti_C0m_0m = CJ0**half
    ti_C1p_1p = CJp
    ti_C1m_1m = CJm 

    aux_list = [ti_C0p_00, ti_C0m_00, ti_C1p_00, ti_C1m_00]
    L2_00 = sum(t**2 for t in aux_list)
    L2_pp = ti_C0p_0p**2 + ti_C1p_1p**2
    L2_mm = ti_C0m_0m**2 + ti_C1m_1m**2
    
    L2_0p = ti_C0p_00*ti_C0p_0p + ti_C1p_00*ti_C1p_1p
    L2_0m = ti_C0m_00*ti_C0m_0m + ti_C1m_00*ti_C1m_1m

    vthep = BOUNDS['vthep']
    vthem = BOUNDS['vthem']

    # prints information to get a good choice of vthep and vthem, depends of VERBOSE and SHOW_ABCDE
    print_abcde(vthep, vthem, L2_00, L2_pp, L2_mm, L2_0p, L2_0m)

    C0 = L2_00 + L2_0p * vthep + L2_0m * vthem
    Cp = L2_pp + L2_0p / vthep
    Cm = L2_mm + L2_0m / vthem

    L_inv_save = BOUNDS['Lstab_inv']
    L_inv_comp = (C0 + max(Cp, Cm))**half
    assert verify.o(L_inv_comp)
    add_check(L_inv_comp < L_inv_save, f"L_inv_stab norm not valid: {L_inv_comp} not less than {L_inv_save}") 

    #######################################
    # --- HARD LEMMAS ---
    Linv_bound  = BOUNDS['Lstab_inv']
    prod_Jfv_u1 = BOUNDS['Jfvu1']
    fvap_L2_sq  = BOUNDS['fvap_stab_L2_sq']
    fvap_H1_sq  = BOUNDS['fvap_stab_H1_sq']
    resi_L2_sq  = BOUNDS['resi_stab_L2_sq']
    Dthe_op     = BOUNDS['Dthe_op']

    rad_stab = BOUNDS['rad_stab']
    ell_stab = BOUNDS['lip_stab']
    eige_lare = BOUNDS['lare']
    
    laap = load_data.eigen_value_approx()

    x_max = prod_Jfv_u1 / CJ0**half
    assert verify.o(x_max)

    fvap_L2 = fvap_L2_sq**half
    fvap_H1 = fvap_H1_sq**half
    resi_L2 = resi_L2_sq**half
    

    aux_nume = Linv_bound * (x_max + fvap_L2)
    aux_a = resi_L2 + Dthe_op * fvap_H1
    aux_b = Dthe_op

    a = aux_a * aux_nume
    b = aux_b * aux_nume
    c = x_max

    qq = ONE
    ll = b - c
    aa = a
    discriminant = ll**2 - 4*qq*aa
    
    
    x_min = - (ll + discriminant**half) * half
    assert verify.o(x_min)

    xx = rad_stab
    q_eta_x = ((xx + fvap_H1) * Dthe_op + resi_L2) / (x_max - xx)
    q_deta_x = ((2*xx +  fvap_H1 + fvap_L2) * Dthe_op + resi_L2) / (x_max - xx)**2

    Tstab0 = Linv_bound * ( (xx + fvap_L2) *  q_eta_x + (xx + fvap_H1) * Dthe_op + resi_L2 )
    Tstab1 = Linv_bound * ( (xx + fvap_L2) * q_deta_x + q_eta_x + Dthe_op)
    eige_real = laap.real() - q_eta_x

    cont_exis = Tstab0 < xx
    lips_cond = Tstab1 < ell_stab
    eige_cond = eige_real > eige_lare

    assert verify.o(Tstab0)
    assert verify.o(Tstab1)
    assert verify.o(q_eta_x)
    
    add_check(discriminant > 0, f"Tstab discriminant not positive: {discriminant}")
    add_check(x_min < x_max, f"x_min: {x_min} is not less than x_max: {x_max}")
    add_check(xx < x_max, f"xx: {xx} is not less than x_max: {x_max}")
    
    add_check(cont_exis, f"Radius for Tstab not valid: {Tstab0} not less than {xx}")           # Lemma 3.12
    add_check(lips_cond, f"Tstab not Lipschitz: {Tstab1}")                                     # Lemma 3.13
    add_check(eige_cond, f"Eigen_real not positive: {eige_real} not greater than {eige_lare}") # Lemma 3.14

    #######################################
    # construct output
    verified = all(ok for ok, _ in checks)
    description = None if verified else "\n\n".join(
        desc for ok, desc in checks if not ok
    )
    return verified, description # `description` is `None` if verified is `True`
    
    


# In[ ]:




