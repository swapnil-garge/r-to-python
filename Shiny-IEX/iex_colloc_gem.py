
#%%
import numpy as np
from scipy.special import sh_jacobi

# ==============================================================================
# Helper Functions for Code Reusability
# ==============================================================================

def _calculate_polynomial_derivatives(roots: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculates the first three derivatives of the Lagrange polynomials at the given roots.
    This helper function abstracts the logic common to both rad_colloc and ax_colloc.

    Args:
        roots: A 1D NumPy array of collocation points.

    Returns:
        A tuple containing three 1D NumPy arrays:
        - p1_derivs: Values related to the first derivative.
        - p2_derivs: Values related to the second derivative.
        - p3_derivs: Values related to the third derivative.
    """
    n_points = len(roots)
    p1_derivs = np.zeros(n_points)
    p2_derivs = np.zeros(n_points)
    p3_derivs = np.zeros(n_points)

    for i in range(n_points):
        x_i = roots[i]
        # Select all other roots
        j_values = np.delete(roots, i)
        delta = x_i - j_values

        # Initialize temporary arrays for the inner loop calculation
        p1 = np.zeros(n_points)
        p2 = np.zeros(n_points)
        p3 = np.zeros(n_points)
        p1[0] = 1.0

        # This loop calculates the derivatives for each root based on all other roots
        for j in range(n_points - 1):
            p1[j+1] = delta[j] * p1[j]
            p2[j+1] = delta[j] * p2[j] + 2 * p1[j]
            p3[j+1] = delta[j] * p3[j] + 3 * p2[j]

        p1_derivs[i] = p1[-1]
        p2_derivs[i] = p2[-1]
        p3_derivs[i] = p3[-1]

    return p1_derivs, p2_derivs, p3_derivs

def _calculate_first_derivative_matrix(roots: np.ndarray, p1_derivs: np.ndarray, p2_derivs: np.ndarray) -> np.ndarray:
    """
    Calculates the first derivative collocation matrix (A).
    This logic is identical for both the radial (Ar) and axial (AZ) cases.

    Args:
        roots: A 1D NumPy array of collocation points.
        p1_derivs: Values related to the first derivative from _calculate_polynomial_derivatives.
        p2_derivs: Values related to the second derivative from _calculate_polynomial_derivatives.

    Returns:
        A 2D NumPy array representing the first derivative matrix.
    """
    n_points = len(roots)
    # Create a matrix of root differences: diff[i, j] = roots[i] - roots[j]
    diff_matrix = roots.reshape(-1, 1) - roots

    # Create a boolean mask to identify the diagonal (where diff_matrix is 0)
    identity_mask = np.eye(n_points, dtype=bool)

    # Calculate off-diagonal elements using vectorized operations
    # This is equivalent to: p1[i] / (p1[j] * (roots[i] - roots[j]))
    A = np.divide(
        p1_derivs.reshape(-1, 1) / p1_derivs,
        diff_matrix,
        where=~identity_mask,
        out=np.zeros((n_points, n_points))
    )

    # Calculate and fill the diagonal elements
    diag_A = 0.5 * p2_derivs / p1_derivs
    np.fill_diagonal(A, diag_A)

    return A

# ==============================================================================
# Main Translated Functions
# ==============================================================================

def rad_colloc(N: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Original R Code:
    rad_colloc <- function(N){
      # ... (full R code as provided in prompt) ...
      return(list(B, W))
    }

    Functional Explanation:
    This function translates the R `rad_colloc` function to Python. It calculates the
    collocation matrix 'B' for the 1-D radial Laplacian and the Gauss-Radau
    quadrature weights 'W' for a given number of collocation points N in a
    symmetric sphere.

    Key Translation Choices:
    - `jacobi.g.recurrences` and `polynomial.roots` from the R `orthopolynom` package
      are replaced by `scipy.special.sh_jacobi`, which directly computes the roots
      of shifted Jacobi polynomials on the interval [0, 1]. This is a more direct
      and standard approach in the Python scientific computing ecosystem.
    - R's `data.frame` is replaced with NumPy arrays, which are more performant and
      idiomatic for numerical computations in Python.
    - Matrix construction loops from R are vectorized using NumPy broadcasting. This
      is significantly more efficient and concise than element-wise loops.
    - The core logic for calculating the derivative values and the final matrices 'B'
      and weights 'W' remains identical to the original R code.
    """
    if N <= 1:
        raise ValueError("Number of collocation points N must be greater than 1.")
        
    # Number of interior collocation points
    N_int = N - 1

    # In R, the roots are calculated from Jacobi polynomial recurrence relations.
    # The Python equivalent finds the roots of the shifted Jacobi polynomial P_n^(2.5, 1.5)
    # directly, which is more efficient. The results are identical.
    # The R code uses rev() on roots; we use np.sort() for a canonical order.
    raw_roots = sh_jacobi(N_int, 2.5, 1.5).roots
    roots = np.concatenate((np.sort(raw_roots), [1.0]))
 
    print(f"Collocation points: {roots}")

    # Calculate polynomial derivative values
    p1_derivs, p2_derivs, p3_derivs = _calculate_polynomial_derivatives(roots)

    # Calculate the first derivative matrix (Ar)
    Ar = _calculate_first_derivative_matrix(roots, p1_derivs, p2_derivs)

    # Calculate the second derivative matrix (Br) using vectorized operations
    identity_mask = np.eye(N, dtype=bool)
    diff_matrix = roots.reshape(-1, 1) - roots
    
    # Off-diagonal elements for Br
    off_diag_term = np.divide(1.0, diff_matrix, where=~identity_mask, out=np.zeros_like(diff_matrix))
    Br = 2 * Ar * (np.diag(Ar).reshape(-1, 1) - off_diag_term)
    
    # Diagonal elements for Br
    diag_Br = (1./3.) * p3_derivs / p1_derivs
    np.fill_diagonal(Br, diag_Br)

    # Calculate the symmetric equivalent matrices, same as in R
    Br_sym = 4 * roots.reshape(-1, 1) * Br + 6 * Ar

    # Calculate quadrature weights (W), same as in R
    a_weight = 2.0
    w_i_prime = 1.0 / (roots * p1_derivs**2)
    W = (1.0 / (a_weight + 1.0)) * w_i_prime / np.sum(w_i_prime)

    # The R function returns a list(B, W). Python returns a tuple.
    return Br_sym, W


def ax_colloc(NZ: int) -> np.ndarray:
    """
    Original R Code:
    ax_colloc <- function(NZ) {
      # ... (full R code as provided in prompt) ...
      return(AZ)
    }

    Functional Explanation:
    This function translates the R `ax_colloc` function. It computes the axial
    discretization matrix 'AZ' using shifted Legendre polynomials for a given
    number of points NZ.

    Key Translation Choices:
    - As with `rad_colloc_py`, `scipy.special.sh_jacobi` is used for finding roots.
      For shifted Legendre polynomials, the parameters are (1.0, 1.0).
    - The R code constructs the roots as `c(0, rev(...), 1)`. The Python equivalent
      is `np.concatenate(([0.0], np.sort(...), [1.0]))`.
    - The core logic for calculating the first derivative matrix 'AZ' is identical
      to the `Ar` calculation in `rad_colloc`, so it has been refactored into the
      `_calculate_first_derivative_matrix` helper function to avoid code duplication.
    """
    if NZ <= 2:
        raise ValueError("Number of collocation points NZ must be greater than 2.")
        
    # Number of interior points

    NZ_int = NZ - 2

    # Get roots of the shifted Legendre polynomial (Jacobi with alpha=1.0, beta=1.0)
    raw_roots = sh_jacobi(NZ_int, 1.0, 1.0).roots
    roots = np.concatenate(([0.0], np.sort(raw_roots), [1.0]))

    # Calculate polynomial derivative values needed for the matrix
    p1_derivs, p2_derivs, _ = _calculate_polynomial_derivatives(roots)

    # Calculate the first derivative matrix (AZ)
    AZ = _calculate_first_derivative_matrix(roots, p1_derivs, p2_derivs)

    return AZ
