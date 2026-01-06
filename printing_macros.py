#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""
This module is used to format prints. 
"""
from sage.rings.real_arb import RealBall
from sage.rings.complex_arb import ComplexBall
from parameters import PRINT_DIGITS 


# In[10]:


def print_RBF(ball):
    # Returns a formatted string "mid ± rad" for RealBall/ComplexBall values.
    if isinstance(ball, (RealBall, ComplexBall)):
        return f"{ball.mid():.{PRINT_DIGITS}g} ± {ball.rad():.{7}g}"
    raise TypeError(f"Expected a RealBall/ComplexBall, got {type(ball)}")


# In[ ]:


def printer_initial(text, len_line=30, symbol="#"):
    # Prints text, separating it from the next print with a line of symbols.
    # Length of the line depends on the text size.
    lon = max(len(text), len_line)
    print(text+"\n" + symbol*lon)


# In[ ]:


def printer(text, len_line=30, symbol="#"):
    # Prints text, separating it from the next print with a line of symbols.
    # Length os fixed.
    lon = len_line
    print(text+"\n" + symbol*lon)


# In[ ]:


def print_QED(text, number_spaces=40):
    # Prints text, use to end the sentence with '☐'.
    lon = max(len(text)+2,number_spaces)
    print(text + "\n" + " "*lon + "☐")


# In[ ]:


greek_alphabet = {
        "alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ",
        "epsilon": "ε", "zeta": "ζ", "eta": "η", "theta": "θ",
        "iota": "ι", "kappa": "κ", "lambda": "λ", "mu": "μ",
        "nu": "ν", "xi": "ξ", "omicron": "ο", "pi": "π",
        "rho": "ρ", "sigma": "σ", "tau": "τ", "upsilon": "υ",
        "phi": "φ", "chi": "χ", "psi": "ψ", "omega": "ω",
        "Alpha": "Α", "Beta": "Β", "Gamma": "Γ", "Delta": "Δ",
        "Epsilon": "Ε", "Zeta": "Ζ", "Eta": "Η", "Theta": "Θ",
        "Iota": "Ι", "Kappa": "Κ", "Lambda": "Λ", "Mu": "Μ",
        "Nu": "Ν", "Xi": "Ξ", "Omicron": "Ο", "Pi": "Π",
        "Rho": "Ρ", "Sigma": "Σ", "Tau": "Τ", "Upsilon": "Υ",
        "Phi": "Φ", "Chi": "Χ", "Psi": "Ψ", "Omega": "Ω",
        "partial": "∂"
    }


# In[ ]:


def print_letter(name, alphabet=greek_alphabet):
    # To use another characters when printing
    return alphabet.get(name, "||CHARACTER NOT FOUND||")

