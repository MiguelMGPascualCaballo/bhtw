#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from sage.all import *
from collections import deque
from parameters import RBF, CBF, ZERO, ONE, ONE_DIV_2, LEN_DOM, LEFT_DOM

fold = "supplementary_data/"

def speed():
    data = _load_tab_data("tw_data.txt")
    return data[0][0]

def coeffs():
    data = _load_tab_data("tw_data.txt")
    return [row[0] for row in data[1:]]

def eigen_value_approx():
    data = _load_tab_data_complex("eigen_data.txt")
    return data[0][0]

def eigen_value_regular():
    data = _load_tab_data_complex("eigen_data.txt")
    return data[1][0]

def polynomial_zeros():
    data = _load_tab_data_complex("polynomial_zeros.txt")
    return [row[0] for row in data]

def theta():
    return 1 / RBF(3)

def eigen_fv():
    """
    Returns Fourier coefficients in order:
    [0,1,...,N-1,-1,...,-N]
    """
    V = ode_stab_sing_stepmat()
    
    len_v1_trunc = V.nrows()-2
    v1_trunc = vector(CBF, len_v1_trunc)

    for kk in range(len_v1_trunc):
        v1_trunc[kk] = V[kk,-1]
    
    return v1_trunc

def Linf_gk_bounds():
    data = _load_tab_data("gk_Linf_bounds.txt")
    return vector(RBF, [x[0] for x in data])

def eigen_u1():
    data = _load_tab_data_complex('sing_u1.txt')
    u1 = vector(CBF, [x[0] for x in data])
    return u1 / u1.norm()

def ode_exis_values():
    return  _ode_values("exis_data.txt")

def ode_stab_sing_values():
    vals_p = _ode_values("sing_p_data.txt")
    vals_m = _ode_values("sing_m_data.txt")
    return  _stab_fixed(vals_p, vals_m) # must return pp,pm,mp,mm,ph,mh

def ode_stab_regu_values():
    vals_p = _ode_values("regu_p_data.txt")
    vals_m = _ode_values("regu_m_data.txt")
    return  _stab_fixed(vals_p, vals_m)

def ode_stab_Jfv_values():
    vals_p = _ode_values("Jfv_p_data.txt")
    vals_m = _ode_values("Jfv_m_data.txt")
    return  _stab_fixed(vals_p, vals_m)
    
def ode_exis_bounds():
    data = _load_tab_data("exis_resis.txt")
    return vector(RBF, [x[0] for x in data])

def ode_stab_sing_bounds(): # TODO
    data = _load_tab_data("sing_resis.txt")
    return vector(RBF, [x[0] for x in data])

def ode_stab_regu_bounds():
    data = _load_tab_data("regu_resis.txt")
    return vector(RBF, [x[0] for x in data])

def ode_stab_Jfv_bounds():
    data = _load_tab_data("Jfv_resis.txt")
    return vector(RBF, [x[0] for x in data])

def ode_exis_stepmat():
    return  matrix(CBF, _load_tab_data_complex('exis_V.txt'))

def ode_stab_sing_stepmat():
    return  matrix(CBF, _load_tab_data_complex('sing_V.txt'))

def ode_stab_regu_stepmat():
    return  matrix(CBF, _load_tab_data_complex('regu_V.txt'))


def load_bounds():
    
    d = {}
    d_str = {}

    with open(fold+'bounds.txt', "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("|", 5)

            try:
                label = parts[0].strip()
                value = parts[2].strip()
                side = parts[4].strip()
            except:
                print(raw)
                print(parts)
                raise Exception("STOP")
            d[label] = RBF(value)
            d_str[label] = [value, label, side]

    return d, d_str

###################################

        
def _load_tab_data(name_file):
    with open(fold+name_file, "r") as f:
        return [list(map(RBF, line.split())) for line in f]

def _load_tab_data_complex(name_file):
    with open(fold+name_file, "r") as f:
        return [[_parse_c_token(tok) for tok in line.split()] for line in f if line.strip()]

def _parse_c_token(tok):
    # "re,im"  -> CBF(RBF(re), RBF(im))
    re_s, im_s = tok.split(",", 1)
    return CBF(RBF(re_s), RBF(im_s))

def _eval_sol(x, values):
    result = CBF(0)
    p = ONE
    for val in values:
        result += val * p
        p *= x
    return result

def _with_leading_one(values):
    """Return homogeneous vector [1, *values] in CBF."""
    v = vector(CBF, len(values) + 1)
    v[0] = ONE
    v[1:] = values
    return v

def _with_leading_zero(values):
    """Return affine vector [0, *values] in CBF."""
    v = vector(CBF, len(values) + 1)
    v[0] = ZERO
    v[1:] = values
    return v

def _new_hom_values_comp(last_val, hom_values, rr):
    hom_values_ = _with_leading_one(hom_values)
    hom_val = _eval_sol(-rr, hom_values_)
    C = last_val / hom_val
    return C * hom_values_
    
def _new_aff_values_comp(last_val, hom_values, aff_values, rr):
    hom_values_ = _with_leading_one(hom_values)
    aff_values_ = _with_leading_zero(aff_values)

    hom_val = _eval_sol(-rr, hom_values_)
    aff_val = _eval_sol(-rr, aff_values_)

    C = (last_val - aff_val) / hom_val
    return C * hom_values_ + aff_values_
    

def interval_from_J(J, len_dom=LEN_DOM, left_dom=LEFT_DOM):
    N, n = J
    len_J = len_dom / (2 ** RBF(N))
    bb = left_dom + len_J * RBF(n)
    aa = bb - len_J
    x0 = (aa + bb) * ONE_DIV_2
    rr = (bb - aa) * ONE_DIV_2
    return aa, bb, x0, rr


# Reading file, for ODE approximations
def _read_exis_blocks(name_file):
    blocks = []
    with open(fold+name_file, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue

            tag, j0, j1 = line.split()
            if tag != "J":
                raise ValueError(f"Expected 'J', got {tag!r}")
            J = (int(j0), int(j1))

            if f.readline().strip() != "VALUES":
                raise ValueError("Expected 'VALUES' after J line")

            values = []
            while True:
                s = f.readline()
                if not s:
                    raise EOFError("Unexpected EOF while reading VALUES block")
                s = s.strip()
                if s == "END":
                    break
                row = [_parse_c_token(t) for t in s.split()]
                values.append(row)

            # optional blank line
            _ = f.readline()

            blocks.append([J, values])
    return blocks

def _ode_values(name_file):
    blocks = _read_exis_blocks(name_file)

    len_coeffs = len(blocks[0][1])

    pos_deque, neg_deque = _split_blocks_by_domain(blocks)

    sols_pos = _const_sols_0_x(
        pos_deque, len_coeffs, sign=+1
    )
    sols_neg = _const_sols_0_x(
        neg_deque, len_coeffs, sign=-1
    )

    return sols_neg + sols_pos

def _stab_fixed(blocks1, blocks2):

    nn = len(blocks1[0][1]) # 2001
    N = (nn-1)//2 # 1000
    
    result = []
    for block1, block2 in zip(blocks1, blocks2):
        J, values1 = block1
        J2, values2 = block2

        assert J == J2, f"Mismatch: {J} vs {J2}"
        
        valpp = values1[:N]
        valpm = values1[N:-1]
        valmp = values2[:N]
        valmm = values2[N:-1]
        
        resorted_values = valpp + valpm + valmp + valmm
        
        resorted_values.append(values1[-1])
        resorted_values.append(values2[-1])
        
        result.append([J, resorted_values])

    return result
    
# Split into pos/neg domains
def _split_blocks_by_domain(blocks):
    """
    Returns (pos_deque, neg_deque) with the same ordering choices you had:
    - pos: appendleft
    - neg: append
    """
    pos = deque()
    neg = deque()

    # you used pop() from a list; we can emulate by iterating reversed
    while blocks:
        blk = blocks.pop()
        J, _ = blk
        _, _, x0, _ = interval_from_J(J)
        if x0.mid() >= 0:
            pos.appendleft(blk)
        else:
            neg.append(blk)

    return pos, neg


# Construct solutions (generic)
def _const_sols_0_x(blocks_deque, len_coeffs, sign):
    """
    sign = +1 for positive domain logic, -1 for negative domain logic.
    This absorbs your +/- rr differences.
    """
    sols = []
    last_hom_val = ONE
    last_aff_values = [ZERO for _ in range(len_coeffs)]

    while blocks_deque:
        J, values = blocks_deque.popleft()
        _, _, _, rr = interval_from_J(J)

        # your code uses rr for pos, -rr for neg in different places
        rr_eff = rr if sign > 0 else -rr

        hom_values = values[0]

        # Homogeneous
        new_hom_values = _new_hom_values_comp(last_hom_val, hom_values, rr_eff)
        last_hom_val = _eval_sol(rr_eff, new_hom_values)

        # Affine pieces; carries the right-endpoint value forward to the next block
        sol_values = [] 
        for idx, aff_values in enumerate(values[1:]):
            last_aff_val = last_aff_values[idx]
            new_aff_values = _new_aff_values_comp(last_aff_val, hom_values, aff_values, rr_eff)
            last_aff_values[idx] = _eval_sol(rr_eff, new_aff_values)
            sol_values.append(new_aff_values)

        sol_values.append(new_hom_values)
        sols.append([J, sol_values])

    return sols





            

