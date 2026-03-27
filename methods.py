#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from sage.all import *
from parameters import RBF, ONE_DIV_2, SQRT_1_DIV_3, ONE_DIV_135, ZERO


# In[ ]:


class integral_1D:
    """
    Callable 1D integrator.

    An instance stores a local integration routine `method` with signature:
        method(x, func) -> RBF
    where `x` is an interval [a, b] (array of RBF) and `func` is a Functions_1D object.
    """

    def __init__(self, method):
        self.method = method
        
        
    @classmethod
    def gauss_lege_2dots(cls):
        """
        2-point Gauss–Legendre enclosure with a remainder term using the 4th derivative:

            ∫_a^b f(t) dt ⊂ rr*(f(x0-rr*w) + f(x0+rr*w)) + rr^5/135 * f^{(4)}([a,b])

        where x0=(a+b)/2, rr=(b-a)/2, w=1/sqrt(3).

        Requirements:
            func.derivatives[4] must exist and accept an interval RBF([a,b]).
        """
        
        def Gauss_Legendre_quad_order2(x, func, half=ONE_DIV_2, w3=SQRT_1_DIV_3, o135=ONE_DIV_135):
            """
            Computes the Gauss–Legendre quadrature enclosure over an interval.

            Variables:
                Input:
                    x:    (RBF) Domain of integration [a, b].
                    func: (Functions_1D) Function to integrate, must provide func.derivatives[4](x).
                Output:
                    result: (RBF) enclosure/approximation of ∫_a^b func(t) dt.
            """
            a, b = x
            X = RBF([a, b])
            x0 = (a + b) * half
            rr = (b - a) * half

            f1 = func(x0 - rr * w3)
            f2 = func(x0 + rr * w3)
            f4 = func.derivatives[4](X)

            return rr * (f1 + f2) + (rr**5) * f4 * o135

        return cls(Gauss_Legendre_quad_order2)

    
    def __call__(self, x, func):
        return self.method(x, func)


# In[ ]:


class image_1D:
    """
    Encloses the image f(x) over an interval x=[a,b] using a Taylor model of order k.

    Stored routine signature:
        method(x, func) -> RBF
    where `func` is usually Functions_1D and provides derivatives up to order k.
    """

    def __init__(self, method, order):
        self.method = method
        self.order = order

    @classmethod
    def taylor_method(cls, kk):
        """
        Taylor enclosure of order kk around the midpoint x0=(a+b)/2.

        For j < kk: uses f^{(j)}(x0) (point evaluation).
        For j = kk: uses f^{(kk)}([a,b]) (interval remainder enclosure).

        Uses power enclosures:
            even j: [0, rr]^j      (since (t-x0)^j >= 0)
            odd  j: [-rr, rr]^j
        """

        def method(x, func, order=kk, half=ONE_DIV_2):
            a, b = x
            X = RBF([a, b])
            x0 = (a + b) * ONE_DIV_2
            rr = (b - a) * ONE_DIV_2

            delta = RBF([-rr, rr])
            delta_pos = RBF([0, rr])

            out = ZERO
            for jj in range(order + 1):
                dj = func.derivatives[order](X) if jj == order else func.derivatives[jj](x0)
                dj /= RBF(factorial(jj))
                dj *= (delta_pos**jj if (jj % 2 == 0) else delta**jj)
                out += dj

            return out
            
        return cls(method, kk)

    def __call__(self, x, func):
        return self.method(x, func)


# In[ ]:




