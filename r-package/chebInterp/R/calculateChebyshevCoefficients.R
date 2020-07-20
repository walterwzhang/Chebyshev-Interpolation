
# calculateChebyshevCoefficients -------------------------------------------------------------------
#' Computes the Chebyshev coefficients from a given function and cheb list
#'
#' Also checks to ensure the cheb$T matrix is orthogonal
#' The rounding down to 0 in the beginning is to account for numerical precision and is controlled
#' by the tolerance parameter
#' The function f only takes one argument
#'
#' @param f Function to be approximated (function)
#' @param cheb List of item from initalizeChebyshevApproximator (list)
#' @param tolerance Numerical Tolerance for rounding down
#' @return A list of Chebyshev coefficients (matrix)
#' @export

calculateChebyshevCoefficients   <-   function(f, cheb, tolerance = 1e-12)
{
    # Assert mutually orthogonal T matrix
    # Round numbers less than the tolerance value to 0 to correct for precision
    diag <- apply(t(cheb$T) %*% (cheb$T),1, function(x) {x[x < tolerance]   <-   0; x})
    stopifnot(all(diag[lower.tri(diag)] == 0, diag[upper.tri(diag)] == 0))

    # Compute Coefficients
    cheb1   <-   cheb
    y       <-   f(cheb1$X)
    if (!any(class(y) == "numeric"))
    {
        y   <-   y[,1]
    }
    cheb1$theta   <-   cheb1$A %*% y
    return(cheb1)
}

# -------------------------------------------------------------------------------------------------
