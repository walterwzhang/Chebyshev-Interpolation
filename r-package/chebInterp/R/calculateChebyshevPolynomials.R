# calculateChebyshevPolynomials -------------------------------------------------------------------
#' Computes the polynomials for a given degree and vector of values.
#'
#' Resultant matrix of polynomials is of size length(x) by N + 1
#'
#' @param x Vector of values to compute the polynomials at (numeric)
#' @param N Highest Degree of the Polynomial (Integer)
#' @return A matrix of the polynomials (matrix)
#' @export

calculateChebyshevPolynomials   <-   function(x, N)
{
    K   <-   length(x)

    # Recursively Compute Polynomials
    T   <-   matrix(0L, ncol = N + 1, nrow = K)
    T[,1]   <-   1
    if (N >= 1)
    {
        T[,2]   <-   x
        if (N >= 2)
        {
            for (k in 2:N)
            {
                T[,k+1]   <-   2*x*T[,k] - T[,k-1]
            }
        }
    }
    return(T)
}

# -------------------------------------------------------------------------------------------------
