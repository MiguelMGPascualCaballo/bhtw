## Data files


### `tw_data.txt`

Plain text file.

- Single row with 1001 space-separated values.
- Value 0: traveling wave speed.
- Values 1–1000: cosine Fourier coefficients of the traveling wave.

This file provides the numerical data defining the approximate traveling wave used throughout the computer-assisted proofs.

### `polynomial_zeros.txt`

Plain text file.

- 1000 rows.
- Each row contains two space-separated values.
- The two values correspond to the real and imaginary parts of a zero of the polynomial $z^N P(z)$, where $N = 1000$.

Only the zeros contained in the open unit disc are listed. These approximations are used as initial guesses for the rigorous root enclosure procedures.




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





