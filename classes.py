#!/usr/bin/env python
# coding: utf-8

# In[1]:


from sage.all import *
from sage.rings.real_arb import RealBall

from printing_macros import print_RBF
from parameters import ZERO, ONE_DIV_2, TWO, RBF, PI


# In[ ]:





# In[2]:


class FourierRealSeries:
    """
    Initialize a real Fourier series with RealBallField (RBF) coefficients.

    Convention:
        - mean is None  <=> exact zero mean
        - mean is not None must be a RealBall (RBF)

    Inputs:
        mean:      RealBall | None
        cos_terms: iterable[RealBall] | None   coefficients [a1, a2, ...] for cos(kx)
        sin_terms: iterable[RealBall] | None   coefficients [b1, b2, ...] for sin(kx)
    """

    def __init__(self, mean=None, cos_terms=None, sin_terms=None):
        """
        Initializes a Fourier series with RBF coefficients.

        Variables:
            Input:
                mean:      (RBF) constant term.
                cos_terms: (list[RBF] | None) cosine coefficients [a_1, a_2, ...].
                sin_terms: (list[RBF] | None) sine coefficients   [b_1, b_2, ...].
            Output:
                self: FourierRealSeries instance with copied coefficient lists.

        Notes:
            - This constructor is strict: all coefficients must already be RealBall (RBF).
              (If you want to accept ints/strings/symbolic and cast to RBF, do it outside
              or replace the isinstance checks with `RBF(val)` casts.)
        """
        
        # Define default coefficient lists if none provided
        cos_terms = cos_terms if cos_terms is not None else []
        sin_terms = sin_terms if sin_terms is not None else []

        #######################################
        # Preallocates variables
        self.cos_terms = [ZERO for _ in range(len(cos_terms))]
        self.sin_terms = [ZERO for _ in range(len(sin_terms))]

        #######################################
        # Defines mean
        # Mean: allow None (exact zero), otherwise require RealBall
        if mean is not None and not isinstance(mean, RealBall):
            raise TypeError("mean must be a RealBall (RBF) or None (exact zero).")
        self.mean = mean

        #######################################
        # Validate and copy cosine coefficients (shallow copy of container)
        cos_list = list(cos_terms)
        for kk, co in enumerate(cos_list, start=1):
            if not isinstance(co, RealBall):
                raise TypeError(f"cos_terms[{kk}] must be a RealBall (RBF).")
        self.cos_terms = cos_list
    
        # Validate and copy sine coefficients (shallow copy of container)
        sin_list = list(sin_terms)
        for kk, si in enumerate(sin_list, start=1):
            if not isinstance(si, RealBall):
                raise TypeError(f"sin_terms[{kk}] must be a RealBall (RBF).")
        self.sin_terms = sin_list

    @classmethod
    def zero(cls, cos_terms=None, sin_terms=None):
        """
        Construct a Fourier series with exact zero mean.

        By convention:
        - mean == None represents exact zero.
        - cos_terms and sin_terms define the non-constant part.
        """
        return cls(
            None,
            cos_terms or [],
            sin_terms or []
        )

    @classmethod
    def constant(cls, c):
        """
        Construct a constant Fourier series.

        Convention:
        - c == None represents the exact zero constant.
        - c != None is treated as a genuine RealBallField element.
        """
        if c is None:
            # Exact zero constant
            return cls(None, [], [])
        return cls(c, [], [])

    @classmethod
    def cosine(cls, coeffs):
        """
        Construct a purely cosine Fourier series
        with exact zero mean.
        """
        return cls(None, list(coeffs), [])

    @classmethod
    def sine(cls, coeffs):
        """
        Construct a purely sine Fourier series
        with exact zero mean.
        """
        return cls(None, [], list(coeffs))

    
    def __add__(self, other):
        """
        Adds Fourier series (or adds a scalar RealBall to the mean term).

        If other is FourierRealSeries:
            (f + g)(x) = (mean_f + mean_g)
                       + Σ_{n>=1} ( (a_n^f + a_n^g) cos(nx) + (b_n^f + b_n^g) sin(nx) ).

        If other is RealBall:
            (f + c)(x) is obtained by shifting the mean term: mean <- mean + c.

        Variables:
            Input:
                other: (FourierRealSeries | RBF)
            Output:
                result: (FourierRealSeries)
        """
        
        # Scalar addition: only allowed for RBF

        if isinstance(other, RealBall):
            return FourierRealSeries.constant(other) + self

        if isinstance(other, FourierRealSeries):
            #######################################
            # Result length: pad to the max sizes
            max_len_c = max(len(self.cos_terms), len(other.cos_terms))
            max_len_s = max(len(self.sin_terms), len(other.sin_terms))
            
            # Preallocate results
            new_mean = FourierRealSeries._add_means(self.mean, other.mean)
            new_cos = [ZERO for _ in range(max_len_c)]
            new_sin = [ZERO for _ in range(max_len_s)]
    
            #######################################
            # Add cosine terms
            for n, co in enumerate(self.cos_terms):
                new_cos[n] += co
            for n, co in enumerate(other.cos_terms):
                new_cos[n] += co
    
            #######################################
            # Add sine terms
            for n, si in enumerate(self.sin_terms):
                new_sin[n] += si
            for n, si in enumerate(other.sin_terms):
                new_sin[n] += si
            
            return FourierRealSeries(new_mean, new_cos, new_sin)

        raise TypeError("Trying to add with not RBF/FourierRealSeries data.")
        
    
    def __radd__(self, other):
        """Adding from the right."""
        return self + other

    
    def __sub__(self, other):
        """
        Subtracts Fourier series:
            f - g = f + (-1)*g

        Variables:
            Input:
                other: (FourierRealSeries)
            Output:
                result: (FourierRealSeries)
        """
        substracter = other * RBF(-1)
        return self + substracter
        
    
    def __mul__(self, other):
        """
        Multiplies two Fourier series (or multiplies by a scalar RealBall).

        Uses trigonometric product-to-sum identities:
            cos(mx)cos(nx) = [cos((m+n)x) + cos((m-n)x)]/2
            sin(mx)sin(nx) = [cos((m-n)x) - cos((m+n)x)]/2
            sin(mx)cos(nx) = [sin((m+n)x) + sin((m-n)x)]/2

        Variables:
            Input:
                other: (FourierRealSeries | RBF)
            Output:
                result: (FourierRealSeries)

        Note: new_cos and new_sin are preallocated to their maximum possible size;
        some trailing coefficients may remain ZERO if no frequency interaction reaches them.
        """

        #######################################
        # Scalar multiplication: allowed only for RealBall
        if isinstance(other, RealBall):
            # If other is a number, transforms it in FourierRealSeries
            return FourierRealSeries.constant(other) * self

        if isinstance(other, FourierRealSeries):
            
            #######################################
            # Compute tight supports
            result_size_c = self._max_cos_mode_with(other)
            result_size_s = self._max_sin_mode_with(other)

            # Mean of the product
            new_mean = FourierRealSeries._mul_means(self.mean, other.mean)

            # Preallocate coefficient arrays
            new_cos = [ZERO for _ in range(result_size_c)]
            new_sin = [ZERO for _ in range(result_size_s)]
            
            #######################################
            # mean term of self contribution
            
            # self.mean * other's Fourier terms
            if self.mean is not None:
                for kk, co in enumerate(other.cos_terms):
                    new_cos[kk] += self.mean * co
                for kk, si in enumerate(other.sin_terms):
                    new_sin[kk] += self.mean * si

            # other.mean * self's Fourier terms
            if other.mean is not None:
                for kk, co in enumerate(self.cos_terms):
                    new_cos[kk] += other.mean * co
                for kk, si in enumerate(self.sin_terms):
                    new_sin[kk] += other.mean * si
            
            #######################################
            # cos(mx)cos(nx) contribution
            for m, a_m in enumerate(self.cos_terms):                    
                for n, a_n in enumerate(other.cos_terms):
                        
                    # cos(mx)cos(nx) = [cos((m+n)x) + cos(|m-n|x)]/2
                    coef = a_m * a_n * ONE_DIV_2
                    
                    # cos((m+n)x) term
                    sum_idx = m + n
                    new_cos[sum_idx+1] += coef
                    
                    # cos(|m-n|)x term
                    diff_idx = abs(m - n)
                    if diff_idx == 0:
                        new_mean = FourierRealSeries._add_means(new_mean, coef)
                    else:
                        new_cos[diff_idx - 1] += coef
            
            #######################################
            # sin(mx)sin(nx) contribution
            for m, b_m in enumerate(self.sin_terms):                    
                for n, b_n in enumerate(other.sin_terms):
                        
                    # sin(mx)sin(nx) = [cos(|m-n|x) - cos((m+n)x)]/2
                    coef = b_m * b_n * ONE_DIV_2
                    
                    # cos(|m-n|)x term
                    diff_idx = abs(m - n)
                    if diff_idx == 0:
                        new_mean = FourierRealSeries._add_means(new_mean, coef)
                    else:
                        new_cos[diff_idx-1] += coef
                    
                    # -cos((m+n)x) term
                    sum_idx = m + n
                    new_cos[sum_idx+1] -= coef
            
            #######################################
            # cos(mx)sin(nx) contribution
            for m, a_m in enumerate(self.cos_terms):                    
                for n, b_n in enumerate(other.sin_terms):
                        
                    # cos(mx)sin(nx) = [sin((m+n)x) + sin((n-m)x)]/2
                    coef = a_m * b_n * ONE_DIV_2
                    
                    # sin((m+n)x) term
                    sum_idx = m + n
                    new_sin[sum_idx+1] += coef
                    
                    # sin((n-m)x) or -sin((m-n)x) term
                    if n > m:
                        diff_idx = n - m
                        new_sin[diff_idx-1] += coef
                    elif m > n:
                        diff_idx = m - n
                        new_sin[diff_idx-1] -= coef
            
            #######################################
            # sin(mx)cos(nx) contribution
            for m, b_m in enumerate(self.sin_terms):                    
                for n, a_n in enumerate(other.cos_terms):
                        
                    # sin(mx)cos(nx) = [sin((m+n)x) + sin((m-n)x)]/2
                    coef = b_m * a_n * ONE_DIV_2
                    
                    # sin((m+n)x) term
                    sum_idx = m + n
                    new_sin[sum_idx+1] += coef
                    
                    # sin((m-n)x) or -sin((n-m)x) term
                    if m > n:
                        diff_idx = m - n
                        new_sin[diff_idx-1] += coef
                    elif n > m:
                        diff_idx = n - m
                        new_sin[diff_idx-1] -= coef

            return FourierRealSeries(new_mean, new_cos, new_sin) 
        
        raise TypeError("Trying to multiply with not RBF/FourierRealSeries data.")
    
    def __rmul__(self, other):
        """Multiplying from the right."""
        return self * other

    @staticmethod
    def _add_means(a, b):
        """Add two means, where None represents exact zero."""
        if a is None:
            return b
        if b is None:
            return a
        return a + b

    @staticmethod
    def _mul_means(a, b):
        """Multiply two means, where None represents exact zero."""
        if a is None or b is None:
            return None
        return a * b

    def _max_cos_mode_with(self, other):
        """
        Maximum k >= 1 such that cos(kx) may appear in the product.
        Returns 0 if no cosine terms can appear.
        """
        mc = len(self.cos_terms)
        ms = len(self.sin_terms)
        nc = len(other.cos_terms)
        ns = len(other.sin_terms)

        self_has_mean  = self.mean is not None
        other_has_mean = other.mean is not None

        candidates = []

        # mean * cos
        if other_has_mean and mc > 0:
            candidates.append(mc)
        if self_has_mean and nc > 0:
            candidates.append(nc)

        # cos * cos
        if mc > 0 and nc > 0:
            candidates.append(mc + nc)

        # sin * sin -> cos
        if ms > 0 and ns > 0:
            candidates.append(ms + ns)

        return max(candidates) if candidates else 0

    def _max_sin_mode_with(self, other):
        """
        Maximum k >= 1 such that sin(kx) may appear in the product.
        Returns 0 if no sine terms can appear.
        """
        mc = len(self.cos_terms)
        ms = len(self.sin_terms)
        nc = len(other.cos_terms)
        ns = len(other.sin_terms)

        self_has_mean  = self.mean is not None
        other_has_mean = other.mean is not None

        candidates = []

        # mean * sin
        if other_has_mean and ms > 0:
            candidates.append(ms)
        if self_has_mean and ns > 0:
            candidates.append(ns)

        # sin * cos, cos * sin
        if ms > 0 and nc > 0:
            candidates.append(ms + nc)
        if mc > 0 and ns > 0:
            candidates.append(mc + ns)

        return max(candidates) if candidates else 0



    
    def __call__(self, x):
        """
        Evaluates the Fourier series at x.

        Variables:
            Input:
                x: (RBF) evaluation point.
            Output:
                result: (RBF) enclosure of f(x).
        """
        result = self.mean if self.mean is not None else ZERO
        for kk, co in enumerate(self.cos_terms, start=1):
            result += co * cos(kk * x)
        for kk, si in enumerate(self.sin_terms, start=1):
            result += si * sin(kk * x)
        return result

    
    def dx(self):
        """
        Computes the derivative of the Fourier series.

        Uses:
            d/dx[a_n*cos(nx)] = -n*a_n*sin(nx)
            d/dx[b_n*sin(nx)] =  n*b_n*cos(nx)

        Variables:
            Output:
                result: (FourierRealSeries) derivative.
        """
        if not self.cos_terms and not self.sin_terms:
            return FourierRealSeries.zero()
        
        # Maximum length
        c_len = len(self.sin_terms)
        s_len = len(self.cos_terms)
        
        # Preallocate arrays
        new_cos = [ZERO for _ in range(c_len)]
        new_sin = [ZERO for _ in range(s_len)]
        
        #######################################
        # Derivative of cos terms -> sin terms
        for n, co in enumerate(self.cos_terms):
            new_sin[n] = -(n+1) * co
        
        #######################################
        # Derivative of sin terms -> cos terms
        for n, si in enumerate(self.sin_terms):
            new_cos[n] = (n+1) * si
        
        return FourierRealSeries(None, new_cos, new_sin)
        

    def ix(self):
        """
        Computes an antiderivative of the Fourier series with zero mean.

        Uses:
            dx^{-1}[a_n*cos(nx)] =  (a_n/n)*sin(nx)
            dx^{-1}[b_n*sin(nx)] = -(b_n/n)*cos(nx)

        Variables:
            Output:
                result: (FourierRealSeries) antiderivative with mean = 0.
        """
        if not self.cos_terms and not self.sin_terms:
            return FourierRealSeries.zero()
        
        # Maximum length
        c_len = len(self.sin_terms)
        s_len = len(self.cos_terms)
        
        # Preallocate arrays
        new_cos = [ZERO for _ in range(c_len)]
        new_sin = [ZERO for _ in range(s_len)]
        
        #######################################
        # Anti-derivative of cos terms -> sin terms
        for n, co in enumerate(self.cos_terms):
            new_sin[n] = co / RBF(n+1)
        
        #######################################
        # Anti-derivative of sin terms -> cos terms
        for n, si in enumerate(self.sin_terms):
            new_cos[n] = -si / RBF(n+1)
        
        return FourierRealSeries(None, new_cos, new_sin)
        

    def hx(self):
        """
        Computes the Hilbert transform of a Fourier series.

        Uses:
            H[a_n*cos(nx)] =  a_n*sin(nx)
            H[b_n*sin(nx)] = -b_n*cos(nx)

        Variables:
            Output:
                result: (FourierRealSeries) Hilbert transform.
        """
        if not self.cos_terms and not self.sin_terms:
            return FourierRealSeries.zero()
        
        # Maximum length
        c_len = len(self.sin_terms)
        s_len = len(self.cos_terms)
        
        # Preallocate arrays
        new_cos = [ZERO for _ in range(c_len)]
        new_sin = [ZERO for _ in range(s_len)]
        
        #######################################
        # Hilbert transform of cos terms -> sin terms
        for n, co in enumerate(self.cos_terms):
            new_sin[n] = co
        
        #######################################
        # Hilbert transform of sin terms -> cos terms
        for n, si in enumerate(self.sin_terms):
            new_cos[n] = -si
        
        return FourierRealSeries(None, new_cos, new_sin)
        
    
    def norm_homogeneous(self, s=0):
        """
        Computes the homogeneous Sobolev s-norm of the represented function.

        Convention used:
            ||f||_s^2 = (1/2) * Σ_{n>=1} (n^{2s}) (a_n^2 + b_n^2)

        Variables:
            Input:
                s: (int | RBF) Sobolev exponent (default 0).
            Output:
                norm: (RBF) norm value.
        """

        exponent = TWO*s;
        
        #######################################
        # Contribution of cos_terms to the norm
        sum_squares_cos = ZERO
        for n, co in enumerate(self.cos_terms):
            sum_squares_cos += co**2 * (n+1)**exponent
            
        #######################################
        # Contribution of sin_terms to the norm
        sum_squares_sin = ZERO
        for n, si in enumerate(self.sin_terms):
            sum_squares_sin += si**2 * (n+1)**exponent

        sum_squares = (sum_squares_cos + sum_squares_sin) * ONE_DIV_2
        return sum_squares**ONE_DIV_2

    
    def __str__(self, tab_len=4, max_terms=5):
        """
        String representation (truncated preview).

        Variables:
            Input:
                tab_len:   (int) indentation spaces.
                max_terms: (int) number of coefficients to show per list.
            Output:
                text: (str) formatted coefficient summary.
        """
    
        tab_space = " " * tab_len

        #######################################
        # Mean term
        result = "Mean:\n" + tab_space + f"{print_RBF(self.mean) if self.mean is not None else 'zero'}"

        #######################################
        # Cosine terms
        if self.cos_terms:
            cos_preview = self.cos_terms[:max_terms]
            cos_str_aux = "\n".join(tab_space + print_RBF(c) for c in cos_preview)
            result += f"\nCosine terms:\n{cos_str_aux}"
            if len(self.cos_terms) > max_terms:
                result += f"\n{tab_space}...and {len(self.cos_terms) - max_terms} more"

        #######################################
        # Sine terms
        if self.sin_terms:
            sin_preview = self.sin_terms[:max_terms]
            sin_str_aux = "\n".join(tab_space + print_RBF(s) for s in sin_preview)
            result += f"\nSine terms:\n{sin_str_aux}"
            if len(self.sin_terms) > max_terms:
                result += f"\n{tab_space}...and {len(self.sin_terms) - max_terms} more"

        return result

    


# In[16]:


class Functions_1D:
    """
    This class is used to define every 1D function `f` we work with.
    The magic method __call__() evaluates the function, and we store derivatives of `f`
    as other callable objects inside `self.derivatives`.

    Design:
        - self.func is the order-0 evaluator: f(x)
        - self.derivatives is a dictionary:
              derivatives[k] is a callable representing f^{(k)}(x)
          (callables can be lambdas, Functions_1D, FourierRealSeries, etc.)

    Variables:
        Attributes:
            func:        (callable) function evaluator, f(x).
            name:        (str) human-readable name.
            derivatives: (dict[int, callable]) dictionary of derivative evaluators.
                         Must contain order 0 and may contain higher orders.
    """

    def __init__(self, func, name=None):
        """
        Initialize a 1D function.

        Variables:
            Input:
                func: (callable) the function evaluator f(x).
                name: (str | None) optional name for the function.
            Output:
                self: Functions_1D object, with derivatives[0] defined.
        """
        self.func = func
        self.name = name if name is not None else "unnamed_function"
        
        #######################################
        # Dictionary of derivatives of the function.
        # Contract: every value stored in this dictionary must be callable,
        # i.e. it can be evaluated as f(x). This includes lambdas and objects
        # implementing __call__ (such as FourierRealSeries instances).
        self.derivatives = {0:func}

    @classmethod
    def identity(cls, tot_ders=2, name=None):
        """
        Creates the identity function id(x) = x.

        Derivatives:
            id'(x) = 1
            id^(k)(x) = 0 for k >= 2

        Variables:
            Input:
                tot_ders: (int) highest derivative order to store.
                name:     (str | None) optional name.
            Output:
                result: (Functions_1D) identity function with derivatives until order tot_ders.
        """
        if not name:
            name = "id"

        def f0(x):
            return x

        result = cls(f0, name=name)

        #######################################
        # First derivative is constant 1
        def one(x):
            return RBF(1)

        # Higher derivatives are 0
        def zero(x):
            return ZERO

        if tot_ders >= 1:
            result.derivatives[1] = one
        for k in range(2, tot_ders + 1):
            result.derivatives[k] = zero

        return result
    
    @classmethod
    def constant(cls, c, tot_ders, name=None):
        """
        Create a constant function f(x) = c.

        All derivatives of positive order are identically 0.

        Variables:
            Input:
                c:        (RBF) constant value.
                tot_ders:  (int) number of derivatives to store (highest order available).
                name:     (str | None) optional function name.
            Output:
                result: (Functions_1D) constant function with derivatives up to tot_ders.
        """
        if not name:
            name = f"const({print_RBF(c)})"

        def f0(x, c=c):
            return c

        result = cls(f0, name=name)

        #######################################
        # Derivatives: f^(k)=0 for k>=1
        def zero(x): return ZERO
            
        for idx in range(tot_ders):
            result.derivatives[idx + 1] = zero
        
        return result

    
    @classmethod
    def from_FourierRealSeries(cls, serie_in, tot_derivatives_in=2, name=None):
        """
        Create a Functions_1D object from a FourierRealSeries (or any callable series-like object).

        The series itself is used as the order-0 callable, and higher derivatives are obtained
        by repeatedly calling `.dx()` on the previous derivative object.

        Variables:
            Input:
                serie_in:            (FourierRealSeries) callable Fourier series object.
                tot_derivatives_in:  (int) highest derivative order to store.
                name:                (str | None) function name.
            Output:
                result: (Functions_1D) function with derivatives[0..tot_derivatives_in].
        """
        if not name:
            name = "unnamed_fourier_real_series"

        #######################################
        # Initialize object with the series evaluator
        result = cls(serie_in, name)

        #######################################
        # Define derivatives: d^k/dx^k serie_in, using the series .dx() method
        for kk in range(1, tot_derivatives_in+1):
            result.derivatives[kk] = result.derivatives[kk - 1].dx()

        return result
        
    
    def __call__(self, x):
        """
        Evaluate function at point x.

        Variables:
            Input:
                x: (RBF) evaluation point.
            Output:
                value: (RBF) f(x).
        """
        return self.func(x)
    

    def dx(self):
        """
        Creates a new Functions_1D object representing the derivative of `self`.

        This method shifts the stored derivative table:
            new.derivatives[0] = old.derivatives[1]
            new.derivatives[1] = old.derivatives[2]
            ...
        so the new object has one fewer available derivative order.

        Variables:
            Output:
                result: (Functions_1D) derivative function.
        """
        # Invariant: derivatives contains all orders 0..N with no gaps.
        tot_derivatives = len(self.derivatives)-1
        if tot_derivatives <= 0:
            raise Exception('Not enough derivatives defined.')

        name = f"d{self.name}"
        result = Functions_1D(self.derivatives[1], name)

        #######################################
        # Shift derivatives down by one order
        for jj in range(1, tot_derivatives):
            result.derivatives[jj] = self.derivatives[jj+1]
            
        return result

    
    def __add__(self, other):
        """
        Adds two functions, or adds a scalar to a function.

        If other is a scalar (RealBall), we convert it into a constant function
        with the same derivative depth, then add.

        If other is Functions_1D, we add order-by-order derivatives.

        Variables:
            Input:
                other: (Functions_1D | RealBall)
            Output:
                result: (Functions_1D) sum function.
        """

        #######################################
        # Scalar case: convert to constant function and reuse logic
        if isinstance(other, RealBall):
            return self + Functions_1D.constant(other, len(self.derivatives)-1)

        #######################################
        # Function + Function case
        if isinstance(other, Functions_1D):
            
            # Creates the function name and function
            name = f"{self.name} + {other.name}"
            func = lambda x: self.func(x) + other.func(x)

            # Initializes the instance
            result = Functions_1D(func, name)
                                  
            #######################################
            # Transfer derivatives (up to the minimum available order)
            max_order = min(len(self.derivatives), len(other.derivatives))-1
            for order in range(1, max_order+1):
                f_deriv = self.derivatives[order]
                g_deriv = other.derivatives[order]
                result.derivatives[order] = lambda x, f=f_deriv, g=g_deriv: f(x) + g(x)

            return result
            
        raise TypeError(
            f"'other' is not a RealBall or Functions_1D. "
            f"Trying to add a {type(other)}."
        )

    
    def __mul__(self, other):
        """
        Multiplies this function with another function or scalar.

        If other is a scalar (RealBall), we convert it into a constant function
        with the same derivative depth, then multiply.

        If other is Functions_1D, we compute derivatives using Leibniz rule:
            (fg)^{(n)} = Σ_{j=0}^n binom(n,j) f^{(j)} g^{(n-j)}.

        Variables:
            Input:
                other: (Functions_1D | RealBall)
            Output:
                result: (Functions_1D) product function.
        """

        #######################################
        # Scalar case: convert to constant function and reuse logic
        if isinstance(other, RealBall):
            return self * Functions_1D.constant(other, len(self.derivatives) - 1)
                
        #######################################
        # Function * Function case
        if isinstance(other, Functions_1D):

            # Creates the function name and function
            name = f"{self.name} * {other.name}"
            func = lambda x: self.func(x) * other.func(x)

            # Initializes the instance
            result = Functions_1D(func, name)
        
            #######################################
            # Computes derivatives using Leibniz rule (generalized product rule)
            fders = self.derivatives
            gders = other.derivatives
            
            max_order = min(len(fders), len(gders)) - 1
            
            # Precompute Leibniz recipes:
            # for each n, store list of (binom(n,j), f^{(j)}, g^{(n-j)})
            recipes = {}
            for n in range(1, max_order + 1):
                terms = []
                for j in range(0, n + 1):
                    terms.append((binomial(n, j), fders[j], gders[n - j]))
                recipes[n] = terms
        
            for kk in range(1, max_order + 1):
                terms = recipes[kk]
        
                def kth_derivative(x, terms=terms):
                    val = ZERO
                    for coeff, f_j, g_nj in terms:
                        val += coeff * f_j(x) * g_nj(x)
                    return val
        
                result.derivatives[kk] = kth_derivative

            return result
                
        raise TypeError(
            f"'other' is not a RealBall or Functions_1D. "
            f"Trying to multiply by a {type(other)}."
        )

    
    def __rmul__(self, other):
        """ Right multiplication (for scalar * function)."""
        return self.__mul__(other)

    
    def __radd__(self, other):
        """ Right add (for scalar + function)."""
        return self.__add__(other)

    
    def __sub__(self, other):
        """
        Subtracts scalars or functions:
            f - g = f + (-1) * g

        Variables:
            Input:
                other: (Functions_1D | RealBall)
            Output:
                result: (Functions_1D) difference function.
        """
        result = self + other * RBF(-1)

        if isinstance(other, RealBall):
            name = f"{self.name} - {print_RBF(other)}"
        elif isinstance(other, Functions_1D):
            name = f"{self.name} - {other.name}"
        else:
            name = f"{self.name} - ({type(other)})"

        result.name = name
        return result

    
    def __pow__(self, n):
        """
        Power by a natural exponent n (0 not allowed).

        Uses exponentiation by squaring:
            f^n = Π f^{(bit)}.

        Variables:
            Input:
                n: (int) natural exponent, n >= 1.
            Output:
                result: (Functions_1D) function representing (self)^n.
        """

        result = Functions_1D.constant(RBF(1), len(self.derivatives)-1, name=f"({self.name})^{n}")
        base = self
            
        while n > 0:
            if n & 1:
                result = result * base
            n >>= 1
            if n:
                base = base * base

        return result


    def inverse(self, name=None):
        """
        Returns 1/self as a Functions_1D object with derivatives up to the available order.

        Derivatives are computed via the Faà di Bruno formula specialized to the inverse:
            (1/f)^{(n)} = sum_{partitions of n} n! / (prod n_j! * (j!)^{n_j})
                          * (-1)^{|ns|} * |ns|! / f^{|ns|+1} * prod f^{(j)}^{n_j}
        where the sum is over tuples (n_1,...) with sum j*n_j = n and |ns| = sum n_j.
        Recipes are cached per order to avoid recomputation across evaluations.
        """
        if not name:
            name = f"({self.name})^-1"
    
        result = Functions_1D(lambda x: 1 / self.func(x), name)

        #######################################
        # Helper: generate tuples (n1, ..., nk) s.t. sum(j * n_j) = k
        def generate_partition_counts(k):
            results = []
            def helper(current, j, remaining):
                if j > k:
                    if remaining == 0:
                        results.append(tuple(current))
                    return
                max_count = remaining // j
                for count in range(max_count + 1):
                    helper(current + [count], j + 1, remaining - j * count)
            helper([], 1, k)
            return results
    
        max_order = len(self.derivatives) - 1
    
        #######################################
        # Cache fully-precomputed recipes for each n:
        # recipe[n] = list of (prefactor, abs_ns, nonzero_pairs)
        inverse_recipes = {} 
    
        def get_inverse_recipe(n):
            rec = inverse_recipes.get(n)
            if rec is not None:
                return rec

            partitions = generate_partition_counts(n)
            out = []
            n_fact = factorial(n)
    
            for ns in partitions:
                # abs_ns = sum_j n_j
                abs_ns = sum(ns)
    
                # FIRST: n! / Π_j (n_j! * (j!)^{n_j})
                pref = RBF(n_fact)
                for j, nj in enumerate(ns, start=1):
                    if nj == 0:
                        continue
                    pref /= RBF(factorial(nj) * (factorial(j) ** nj))
    
                # SECOND constant part: (-1)^{abs_ns} * abs_ns!
                pref *= RBF(((-1) ** abs_ns) * factorial(abs_ns))
    
                # THIRD: only keep j where n_j > 0
                nonzero = [(j, nj) for j, nj in enumerate(ns, start=1) if nj > 0]
    
                out.append((pref, abs_ns, nonzero))
    
            inverse_recipes[n] = out
            return out

        #######################################
        # Define derivatives of 1/f using cached recipes
        for kk in range(1, max_order + 1):
            recipe = get_inverse_recipe(kk)
    
            def kth_derivative(x, n=kk, recipe=recipe):
                fx = self(x)  # compute once
                val = ZERO
    
                for pref, abs_ns, nonzero in recipe:
                    term = pref / (fx ** (abs_ns + 1))
                    for j, nj in nonzero:
                        term *= (self.derivatives[j](x) ** nj)
                    val += term
    
                return val
    
            result.derivatives[kk] = kth_derivative
    
        return result

    
    def __str__(self):
        """String representation of the function."""
        return f"Function: {self.name}"

