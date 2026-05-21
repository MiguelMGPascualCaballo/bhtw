# Computer-assisted verification for *Linear instability of a Burgers--Hilbert traveling wave*

This repository accompanies the paper [*Linear instability of a Burgers--Hilbert traveling wave*](https://arxiv.org/abs/2605.03920) and contains the code and supplementary data used to reproduce the computer-assisted parts of the proof.

## Requirements

The code must be run in a SageMath-enabled Python environment. In particular, the following import must work:

```python
from sage.all import *
```

## Note

Some large data files listed below are not included in the repository and must be downloaded from Zenodo:

<https://doi.org/10.5281/zenodo.19250315>

These data files must be placed in `supplementary_data/`.

## Python source files

### `script.py`

Main script reproducing all computer-assisted results.

- Calls the verification routines defined in `lemmas.py`.
- Outputs timing information and verification status for each lemma and the final routine.

Running this script reproduces all computer-assisted checks reported in the paper.

### `lemmas.py`

Each function in this file corresponds to a specific lemma, except for the final routine, which verifies that all substitutions made in the paper are correct.

### `auxiliar_funcs.py`

Generic and reusable adaptive verification algorithms.

### `methods.py`

Reusable numerical methods independent of any specific lemma.

### `explicit_funcs.py`

Functions and routines specific to the problem studied in the paper.

### `verify.py`

Helper functions used to verify that everything works correctly and avoid silent bugs.

### `parameters.py`

Defines all global numerical parameters used in the computations.

- `VERBOSE` controls the amount of printed information.

### `printing_macros.py`

Helper functions to print progress and status of the proof.

### `load_data.py`

Routines for reading the data files in `supplementary_data/`.

### `classes.py`

Utility classes used to simplify parts of the code.

- `FourierRealSeries`
- `Functions_1D`

## Data files

All data files are located in the `supplementary_data/` directory. They are plain text files read by the routines in `load_data.py`.

### `tw_data.txt`

- Read with `_load_tab_data()`.
- The first value is the traveling wave speed.
- The remaining values are the traveling wave coefficients.

### `*_resis.txt`

(e.g. `exis_resis.txt`, `sing_resis.txt`, `regu_resis.txt`, `Jfv_resis.txt`)

- Read with `_load_tab_data()`.
- Upper error bounds for the difference between certain functions of the paper and their numerical approximations.

### `gk_Linf_bounds.txt`

- Read with `_load_tab_data()`.

### `eigen_data.txt`

- Read with `_load_tab_data_complex()`.
- The first line contains \( \lambda^{\rm ap} \).
- The second line contains \( 1.035i \).

### `*_V.txt`

(e.g. `exis_V.txt`, `sing_V.txt`, `regu_V.txt`)

- Read with `_load_tab_data_complex()`.
- Each file contains an approximation of the \( V \) factor in the singular value decomposition of a corresponding matrix.
- The entries are provided with up to \( 30 \) decimal digits of precision.
- `sing_V.txt` is also used in the construction of the function \( f^{\rm ap} \).
- Not included in the repository, available on Zenodo.

### `sing_u1.txt`

- Read with `_load_tab_data_complex()`.
- Approximation of the left singular vector \( u_1 \) associated to the smallest singular value of \( M^{\rm tor}_{\rm stab} \). 

### `polynomial_zeros.txt`

- Read with `_load_tab_data_complex()`.
- Approximation of the zeros of \( z^N P(z) \) contained in the open unit disc.

### `exis_data.txt`

- Read with `_ode_values()`.
- Not included in the repository, available on Zenodo.

### `*_p_data.txt` and `*_m_data.txt`

(e.g. `sing_p_data.txt`, `sing_m_data.txt`, `regu_p_data.txt`, `regu_m_data.txt`, `Jfv_p_data.txt`, `Jfv_m_data.txt`)

- Read with `_ode_values()` and restructured with `_stab_fixed()`.
- Not included in the repository, available on Zenodo.

### `bounds.txt`

- Read with `load_bounds()`.
- Contains all remaining bounds used in the paper.

