#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from sage.all import *
from auxiliar_funcs import *
from classes import *
from explicit_funcs import *
from methods import image_1D, integral_1D
from printing_macros import *
from parameters import *


# In[ ]:


def lemma_norm_xi_L2(speed, coeffs, bound_in):
    """
    This function encloses the L2 norm of the approximated traveling wave, in the text we call it xi:
        xi = c v' + Hv + v'v.

    The function verifies whether the computed norm is bounded by `bound_in`.

    Computations are done symbolically.

    Variables:
        Input:
            speed:    (RBF) speed of the traveling wave: c.
            coeffs:   (RBF array) cosine Fourier coefficients of the traveling wave, v = sum 'coeffs[k]' cos(k x).
            bound_in: (str convertible to RBF) lemma bound to verify against.
        Output:
            verified:   (bool) True if the statement ||xi||_L2 < 'bound_in' is true, otherwise False.
            norm_xi_L2: (RBF) L2 norm of the residue xi.     
    """
    
    # Verifies that speed and coeffs are RBF
    verify.INPUT(speed, coeffs)

    # Transforms str to RBF
    bound = RBF(bound_in)

    #######################################
    # Computes xi
    xi = equ_tw_eval_symbolic(speed, coeffs)
    
    # The norm we use (H0odd) for a real odd function is based on the positive complex modes:
    #     ||xi|| = ( sum_{k>=1} |xi_k|^2 )^(1/2),
    # where xi_k are the complex Fourier coefficients in xi(x) = sum_{k∈Z} xi_k e^{ikx}.
    #
    # For a real odd xi, coefficients satisfy xi_{-k} = -conj(xi_k), hence
    #     sum_{k∈Z} |xi_k|^2 = 2 * sum_{k>=1} |xi_k|^2.
    #
    # `xi.norm_homogeneous()` returns ( sum_{k∈Z} |xi_k|^2 )^(1/2),
    # so to obtain ( sum_{k>=1} |xi_k|^2 )^(1/2) we divide by sqrt(2).
    norm_xi_L2 = xi.norm_homogeneous() / (RBF(2) ** (RBF(1)/RBF(2)))

    #######################################
    # Verifies that no precision was lost during the computation
    verify.OUTPUT(norm_xi_L2)
    
    verified = norm_xi_L2 < bound
    print(f"Statement is {verified}.")

    return verified, norm_xi_L2


# In[ ]:


def lemma_enclose_roots(speed, coeffs, roots_in, radii_roots, verbose=VERBOSE):
    """
    Verifies a root-enclosure lemma for the polynomial associated with the traveling wave,
        P(z) = z^N (`speed` + sum_{j=1}^N (z^j + z^{-j}) * coeffs_j/2).

    The function:
      1) Builds the polynomial P associated to (speed, coeffs).
      2) For each approximate root in `roots_in`, proves that there exists a true root of P
         inside the closed disc D(root, radii_roots).
      3) Checks that all such discs are pairwise disjoint.
      4) Checks that all such discs are contained in the open unit disc.

    Computations use ball arithmetic (RBF/CBF) and precomputed coefficient-derivative bounds.

    Inputs:
        speed:       (RBF) wave speed c.
        coeffs:      (RBF array) cosine Fourier coefficients of v.
        roots_in:    (list of CBF) approximate roots of P.
        radii_roots: (RBF) common radius r used to thicken each root into a disc.
        verbose:     (int) print progress messages.

    Output:
        (bool) True iff all of the following hold:
            (A) Every disc D(root, r) is proven to contain a root of P.
            (B) Every disc D(root, r) is contained in the unit disc (|z|<1).
            (C) The discs D(root, r) are pairwise disjoint.
    """
    
    # Verifies that speed and coeffs are RBF
    verify.INPUT(speed, coeffs)

    # Work on a copy, since we destructively pop roots during the proof loop
    roots=roots_in.copy()
    
    #######################################
    # Construct the polynomial P associated to the traveling wave
    rk = len(coeffs) # N
    polynom = palipoly_coeffs_constructor(speed, coeffs)

    #######################################
    # Precompute auxiliary bounds for derivatives of P on |z|<=1.
    #
    # We build numbers `big_bounds[d]` such that, for each derivative order d,
    #     sup_{|z|<=1} |P^{(d)}(z)| <= big_bounds[d].
    #
    # For |z|<=1 and P(z)=sum_k a_k z^k, we have
    #     |P^{(d)}(z)| <= sum_{k>=d} |a_k| * k*(k-1)*...*(k-d+1),
    # hence the coefficient-based bound below.

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

    #######################################
    # Proof loop:
    # - condA: every proposed disc contains a root (validated via `root_polynomial_is_enclosed`)
    # - circ_dist: minimal distance margin to the unit circle (positive => all discs inside |z|<1)
    # - pair_dist: minimal separation margin between discs (positive => discs are disjoint)
    
    circ_dist = RBF(1)
    pair_dist = RBF(2)

    condA = True
    while roots and condA:
        # Take one approximate root and thicken it into a closed disc of radius `radii_roots`
        root = roots.pop()
        root_enclo = root.add_error(radii_roots)

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
        for other_root in roots:
            other_root_enclo=other_root.add_error(radii_roots)
            len_pair = (root_enclo-other_root_enclo).abs()
            pair_dist = pair_dist.min(len_pair)
            

        if verbose and len(roots) % 50 == 49:
            print(len(roots), len_dis, circ_dist)

    # Final conditions:
    # - circ_dist > 0  => every disc is strictly inside the unit disc
    # - pair_dist > 0  => discs are pairwise disjoint
    condB = circ_dist > ZERO
    condC = pair_dist > ZERO

    verified = condA and condB and condC
    print(f"statement is {verified}.")

    return verified
       


# In[ ]:


def lemma_roots_compatible_exis(speed, coeffs, roots_in, radii_roots):

    """
    Verifies a nondegeneracy condition for a list of enclosed roots.

    Given:
      - the polynomial P associated to (speed, coeffs),
      - approximate roots `roots_in`,
      - and a common enclosure radius `radii_roots`,

    the function thickens each root into a complex ball (disc) and checks that, for every such disc,
    a certain interval-valued expression cannot be an integer number.

    More precisely, letting N = number of roots and z be an enclosed root, we form
        image_root(z) = - z^(N-1) / P'(z),.

    The lemma check is:
      - it is OK if image_root is *not* exactly real, OR
      - it is OK if (even though it might be real) its real part cannot be an integer.

    If for some enclosed root both of the following could happen simultaneously:
      - Im(image_root) contains 0 exactly, and
      - Re(image_root) contains an integer,
    then the lemma fails and the function returns False.

    Inputs:
        speed:       (RBF) wave speed c.
        coeffs:      (RBF array) coefficients defining the traveling wave (used to build P).
        roots_in:    (list of CBF) approximate roots of P.
        radii_roots: (RBF) radius used to enclose each root into a disc.

    Output:
        verify: (bool) True if the compatibility condition holds for all enclosed roots, otherwise False.
    """

    # Verifies that speed and coeffs are RBF
    verify.INPUT(speed, coeffs)
    
    # Work on a copy (avoid mutating the user-provided list)
    roots=roots_in.copy()
    
    #######################################
    # Construct the polynomial P associated to the traveling wave
    polynom = palipoly_coeffs_constructor(speed,coeffs)
   
    #######################################
    # Enclose each approximate root into a complex ball of radius `radii_roots`
    enclosed_roots = []
    for root in roots:
        enclosed_root = root.add_error(radii_roots)
        enclosed_roots.append(enclosed_root)

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
        image_root = - root**(len(roots)-1) / polynom_eval(root,polynom,idx_der=1)

        # Separate real and imaginary parts (both are interval/ball quantities)
        image_root_real=image_root.real()
        image_root_imag=image_root.imag()

        # Pass condition:
        #  - either Im(image_root) is provably nonzero (i.e., does NOT contain 0 exactly), OR
        #  - Re(image_root) is provably not an integer (i.e., does NOT contain any integer).
        #
        # Failure happens only if it could be exactly real AND could hit an integer.
        verified = not image_root_imag.contains_exact(0) or not image_root_real.contains_integer()

    print(f"Statement is {verified}.")
    return verified


# In[ ]:


def lemma_roots_compatible_stab(speed, coeffs, roots_in, radii_roots, lamb):

    """
    Verifies a nondegeneracy condition for a list of enclosed roots.

    Given:
      - a polynomial P associated to (speed, coeffs),
      - approximate roots `roots_in`,
      - a common enclosure radius `radii_roots`,
      - and a parameter `lamb`,

    the function thickens each root into a complex ball (disc) and checks that two
    interval-valued expressions cannot be integer real numbers.

    Let N = number of roots and z be an enclosed root. Define
        aux_term(z) = z^(N-1) / P'(z),
    where P'(z) is the first derivative of P evaluated at z (with enclosure).

    Then we form the two "images"
        image_root1(z) = ( lamb - 1) * aux_term(z),
        image_root2(z) = (-lamb - 1) * aux_term(z).

    For each enclosed root z, the lemma requires BOTH images to avoid being an integer real number
    under enclosure arithmetic:
      - it is OK if image_root is *not* exactly real, OR
      - it is OK if (even though it might be real) its real part cannot contain an integer.

    Failure for an image occurs only if it could be exactly real AND could hit an integer, i.e.
      Im(image_root) contains 0 exactly  AND  Re(image_root) contains an integer.

    Inputs:
        speed:       (RBF) wave speed c.
        coeffs:      (RBF array) coefficients defining the traveling wave (used to build P).
        roots_in:    (list of CBF) approximate roots of P.
        radii_roots: (RBF) radius used to enclose each root into a disc.
        lamb:        (RBF/CBF scalar) parameter λ appearing in the stability condition.

    Output:
        verify: (bool) True if the stability compatibility condition holds for all enclosed roots and for both images, otherwise False.
    """

    # Verifies that speed and coeffs are RBF
    verify.INPUT(speed, coeffs)
    
    # Work on a copy (avoid mutating the user-provided list)
    roots=roots_in.copy()
    
    #######################################
    # Construct the polynomial P associated to the traveling wave
    polynom = palipoly_coeffs_constructor(speed,coeffs)

    #######################################
    # Enclose each approximate root into a complex ball of radius `radii_roots`
    enclosed_roots = []
    for root in roots:
        enclosed_root = root.add_error(radii_roots)
        enclosed_roots.append(enclosed_root)

    #######################################
    # Lemma check over all enclosed roots
    #
    # For each enclosed root z we compute aux_term(z) = z^(N-1)/P'(z) and test the two
    # scaled images (λ-1)*aux_term and (-λ-1)*aux_term against the "integer real" obstruction.
    verified = True 
    while enclosed_roots and verified:
        root = enclosed_roots.pop()
        
        # Compute aux_term(z) = z^(N-1) / P'(z) with enclosure arithmetic
        aux_term = root**(len(roots)-1) / polynom_eval(root,polynom,idx_der=1)

        # Two images: (±λ - 1) * aux_term, for each sign of lamb.
        image_root1 = ( lamb - 1) * aux_term
        image_root2 = (-lamb - 1) * aux_term

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
        veri1 = (not image_root_imag1.contains_exact(0)) or (not image_root_real1.contains_integer())
        veri2 = (not image_root_imag2.contains_exact(0)) or (not image_root_real2.contains_integer())

        # We require both images to satisfy the condition for every root
        verified = veri1 and veri2

    print(f"Statement is {verified}.")
    return verified


# In[ ]:


def lemma_beta_maximum(speed, coeffs, bound_in, order=2):
    """
    Verifies a supremum bound for
        beta(x) := 1 / (c + v(x)),
    where v is an approximation of the traveling wave with speed c.

    The function checks whether
        sup_{x in [0, pi]} |beta(x)| < bound_in,
    using an adaptive 1D enclosure algorithm. To improve accuracy, the method
    may use derivatives of beta up to the prescribed order.

    Inputs:
        speed:    (RBF) wave speed c.
        coeffs:   (RBF array) cosine Fourier coefficients of the traveling wave approximation v.
        bound_in: (str convertible to RBF) lemma bound to verify against.
        order:    (int) order of the supremum method (uses derivatives up to this order).

    Outputs:
        verified: (bool) True if sup_{x in [0, pi]} |beta(x)| < bound_in, otherwise False.
    """
    
    # Verifies that speed, coeffs are RBF.
    verify.INPUT(speed, coeffs)

    # Convert the provided bound to RBF
    guess_bound = RBF(bound_in)

    #######################################
    # Construct beta(x) = 1 / (c + v(x)) together with its derivatives up to the prescribed order.
    w_func = w_function_constructor(coeffs, total_derivatives=order)
    beta_func = beta_function_constructor(speed, w_func)

    #######################################
    # Domain construction
    domain = [ZERO, PI]

    # Choose a kth-order enclosure method for the supremum computation
    method = image_1D.kth_order_method_constructor(order)

    #######################################
    # Verify the supremum bound using an adaptive enclosure algorithm
    verified = sup_bound_prev_esti_adapt_1D(domain, beta_func, guess_bound, method)

    print(f"Statement is {verified}.")
    
    return verified


# In[ ]:


def lemma_betamod_maximum(speed, coeffs, bound_in, order=2):
    """
    Verifies a supremum bound for
        betamod(x) := (pi - x) * beta(x)^2,
    where
        beta(x) := 1 / (c + v(x)),
    and v is an approximation of the traveling wave with speed c.

    The function checks whether
        sup_{x in [0, pi]} |betamod(x)| < bound_in
    using an adaptive 1D enclosure algorithm. To improve accuracy, the method
    may use derivatives of beta up to the prescribed order.

    Inputs:
        speed:  (RBF) wave speed c.
        coeffs: (RBF array) cosine Fourier coefficients of the traveling wave approximation v.
        bound_in: (str convertible to RBF) lemma bound to verify against.
        order:  (int) order of the supremum method (uses derivatives up to this order).

    Outputs:
        verified: (bool) True if sup_{x in [0, pi]} |betamod(x)| < bound_in, otherwise False.
    """
    
    # Verifies that speed, coeffs are RBF.
    verify.INPUT(speed, coeffs)

    # Convert the provided bound to RBF
    guess_bound = RBF(bound_in)

    #######################################
    # Construct beta(x) = 1 / (c + v(x)) together with its derivatives up to the prescribed order.
    w_func = w_function_constructor(coeffs, total_derivatives=order)
    beta_func = beta_function_constructor(speed, w_func)

    # Construct betamod(x) = (pi - x) * beta(x)^2 (and derivatives if needed by the method)
    betamod_func = betamod_constructor(beta_func,order)
    
    #######################################
    # Domain construction
    domain = [ZERO, PI]

    # Choose a kth-order enclosure method for the supremum computation
    method = image_1D.kth_order_method_constructor(order)
    
    #######################################
    # Verify the supremum bound using an adaptive enclosure algorithm
    verified = sup_bound_prev_esti_adapt_1D(domain, betamod_func, guess_bound, method)

    print(f"Statement is {verified}.")

    return verified


# In[ ]:


def lemma_betamod2_L2sq(speed, coeffs, bound_in, abs_tol, rel_tol):
    """
    Verifies an L2-squared bound for the beta-derived function
        betamod2(x) := beta'(x) + i * beta(x)^2,
    where
        beta(x) := 1 / (c + v(x)),
    and v is an approximation of the traveling wave with speed c.

    Since beta takes real values, we have
        |betamod2(x)|^2 = |beta'(x) + i*beta(x)^2|^2 = beta'(x)^2 + beta(x)^4.

    The function checks whether
        ∫_{-pi}^{pi} ( beta'(x)^2 + beta(x)^4 ) dx <= bound_in,
    using an adaptive 1D enclosure algorithm for the integral (Gauss–Legendre quadrature).

    Inputs:
        speed:    (RBF) wave speed c.
        coeffs:   (RBF array) cosine Fourier coefficients of the traveling wave approximation v.
        bound_in: (str convertible to RBF) lemma bound to verify against (for the squared L2 norm).

    Outputs:
        verified: (bool) True if the bound holds, otherwise False.
        norm_betamod2_L2_sq: (RBF) enclosure of ∫_{-pi}^{pi} |beta'(x) + i*beta(x)^2|^2 dx.
    """

    # Verifies that speed and coeffs are RBF
    verify.INPUT(speed, coeffs)

    # Convert the provided bound to RBF
    guess_bound = RBF(bound_in)
    abs_tol = RBF(abs_tol)
    rel_tol = RBF(rel_tol)
    
    #######################################
    # Construct beta(x) = 1 / (c + v(x)) together with its derivatives.
    w_func = w_function_constructor(coeffs, total_derivatives=4)
    beta_func = beta_function_constructor(speed, w_func)

    # Nonnegative integrand: |beta'(x) + i*beta(x)^2|^2 = beta'(x)^2 + beta(x)^4
    integrand = betamod2_constructor(speed, coeffs)

    #######################################
    # Construct domain
    domain = [ZERO, PI]

    # Enclosure method for the integral (adaptive Gauss–Legendre)
    method = integral_1D.Gauss_Legendre_quad_order2_constructor()

    # Compute a rigorous enclosure of the integral
    norm_betamod2_L2_sq = 2*integral_adaptive_computation_1D(domain, integrand, method, abs_tol, rel_tol)

    #######################################
    # Verifies that no precision was lost during the computation
    verify.OUTPUT(norm_betamod2_L2_sq)
    
    # Verify the lemma bound (choose <= if that is your statement)
    verified = norm_betamod2_L2_sq <= guess_bound
    
    print(f"Statement is {verified}.")
    return verified, norm_betamod2_L2_sq


# In[ ]:


def lemma_betaiota(speed, coeffs, bound_in, abs_tol):
    """
    Rigorously verifies an L¹([0, π]) bound for the beta–iota expression
        (iota(π) − iota(x)) · beta(x)²,
    using adaptive interval arithmetic.

    The function iota is defined as a primitive of
        diota(x) := (beta'(x))² + beta(x)⁴ ≥ 0,
    i.e.
        iota(x) = ∫₀ˣ diota(y) dy.

    The function beta is defined by
        beta(x) := 1 / (c + v(x)),
    where v is a traveling-wave approximation with wave speed c = `speed`,
    represented via its cosine Fourier coefficients `coeffs`.

    In the proof of the lemma, iota is chosen as an odd primitive (since diota
    is even), which implies:
      - iota is monotone increasing on [0, π],
      - by symmetry,
            || (iota(π) − iota(x)) · beta(x)² ||_{L¹([0, π])}
          = || (iota(x) − iota(−π)) · beta(x)² ||_{L¹([−π, 0])}.

    This implementation verifies the inequality
        || (iota(π) − iota) · beta² ||_{L¹([0, π])} ≤ bound,
    by computing a rigorous upper bound via *adaptive integration*.
    No fixed spatial partition is used.

    Parameters:
        speed (RBF):
            Wave speed c.
        coeffs (RBF array):
            Cosine Fourier coefficients defining the traveling-wave profile v(x).
        bound_in (RBF or convertible):
            Lemma bound to verify.

    Returns:
        verified (bool):
            True if the L¹ inequality is rigorously verified, otherwise False.
        betaiota_L1 (RBF):
            A rigorous enclosure of the computed L¹ quantity.
            
    """

    # Verifies that speed and coeffs are RBF
    verify.INPUT(speed, coeffs)

    # Convert the provided bound to RBF
    guess_bound = RBF(bound_in)

    abs_tol = RBF(abs_tol)

    #######################################
    # Construct beta(x) = 1 / (c + v(x)) together with its derivatives.
    w_func = w_function_constructor(coeffs, total_derivatives=5)
    beta_func = beta_function_constructor(speed, w_func)

    beta_sq = beta_func**2

    # Nonnegative diota := |beta'(x) + i*beta(x)^2|^2 = beta'(x)^2 + beta(x)^4
    diota = betamod2_constructor(speed, coeffs)

    # iota(x) = int_0^x diota(y) dy
    # betaiota_L1 = int_0^pi beta(y)^2 (iota(pi) - iota(y)) dy
    betaiota_L1 = iota_L1_compute_adaptive(diota, beta_sq, abs_tol)

    # Non precision lost
    verify.OUTPUT(betaiota_L1)

    verified = betaiota_L1 <= guess_bound

    print(f"Statement is {verified}.")
    return verified, betaiota_L1

