## Note

The descriptions in this README correspond to a forthcoming version of the project and may not fully match the current contents of the repository. In addition, the documentation is not yet complete and will be expanded in future updates.

## Data files

All data files are located in the `supplementary_data/` directory and are provided as plain text files and readen by `load_data.py` routines.

### `tw_data.txt`

- Read with `_load_tab_data()`.
- The first value is the traveling wave speed.
- The remaining values are the traveling wave coefficients.

### `*_resis.txt`

- Read with `_load_tab_data()`.
- Upper error bounds for the difference between certain functions of the paper and their numerical approximations.

### `gk_Linf_bounds.txt`

- Read with `_load_tab_data()`.

### `eigen_data.txt`

- Read with `_load_tab_data_complex()`.
- The first line contains \( \lambda^{\rm ap} \).
- The second line contains \( 1.035i \).

### `*_V.txt`

- Read with `_load_tab_data_complex()`.
- Each file contains an approximation of the \( V \) factor in the singular value decomposition of a corresponding matrix.
- The entries are provided with up to \( 30 \) decimal digits of precision.
- `sing_V.txt` is also used in the construction of the function \( f^{\rm ap} \).

### `sing_u1.txt`

- Read with `_load_tab_data_complex()`.
- Approximation of the left singular vector \( u_1 \) associated to the smallest singular value of \( M^{\rm tor}_{\rm stab} \). 

### `polynomial_zeros.txt`

- Read with `_load_tab_data_complex()`.
- Approximation of the zeros of \( z^N P(z) \) contained in the open unit disc.

### `exis_data.txt`

- Read with `_ode_values()`.

### `*_p_data.txt` and `*_m_data.txt`

- Read with `_ode_values()` and restructured with `_stab_fixed()`.

### `BOUNDS.txt`

- Read with `load_bounds()`.

## Python source files

### `script.py`

Main execution script reproducing all computer-assisted results.

- Loads the traveling-wave data from `tw_data.txt` and polynomial root approximations from `polynomial_zeros.txt`.
- Calls the verification routines defined in `lemmas.py`.
- Outputs timing information and verification status for each lemma.

Running this script reproduces all computer-assisted checks reported in the paper.
**Important:** this script must be run in a **SageMath-enabled Python environment**, since it imports `sage.all` and uses Sage’s Arb-based ball arithmetic (`RealBallField`, `ComplexBallField`).

### `script.ipynb`

Jupyter notebook version of `script.py`.

- Contains the same workflow as `script.py`, but organized into notebook cells.
- Useful for interactive inspection of intermediate quantities (e.g. partial bounds, progress prints, timings).
- Reproduces the same certified checks by calling the same routines in lemmas.py.

**Important**: the notebook must be executed with a **SageMath kernel** (or a Jupyter setup configured to use Sage), because it relies on `sage.all` and Arb ball arithmetic types.

### `lemmas.py`

High-level computer-assisted lemmas corresponding directly to statements in the paper.

- Verification of norm bounds, supremum bounds, and integral estimates.
- Root enclosure and compatibility checks for the polynomial $z^N P(z)$.
- Wrapper routines that compare rigorous enclosures with the bounds stated in the lemmas.

Each function in this file corresponds to a specific lemma in the Appendix.


### `auxiliar_funcs.py`

Generic and reusable adaptive verification algorithms.

- Adaptive subdivision for supremum bounds.
- Adaptive integration drivers for one-dimensional integrals.
- Input/output sanity checks to prevent loss of rigor or precision.

These routines implement the core logic shared by many computer-assisted lemmas.


### `methods.py`

Local enclosure methods used by adaptive algorithms.

- Taylor model image enclosures of arbitrary order.
- Gauss–Legendre quadrature with rigorous remainder terms.
- Callable wrappers used by adaptive integration and supremum verification routines.

This file contains reusable numerical methods independent of any specific lemma.


### `explicit_funcs.py`

Defines explicit mathematical objects used in the computer-assisted proofs.

- Construction of functions represented as `Functions_1D objects`.
- Explicit adaptive integration routines: `iota_L1_compute_adaptive`.

These routines build the functions and integrands that are later verified by higher-level lemmas.



### `parameters.py`

Defines all global numerical parameters used in the computations.

- Working precision and real/complex ball fields.
- Default absolute and relative tolerances for adaptive algorithms.
- Maximum iteration and subdivision limits.
- Numerical constants used throughout the code.
  
This file centralizes all numerical settings to ensure consistency and reproducibility.


### `printing_macros.py`

Utility functions for formatted output.

- Orientative printing of real and complex ball numbers.
- Helper functions for displaying intermediate and final results.
- This file has no mathematical content and is only used for readable output.
  
This file has no mathematical content and is only used for readable output





