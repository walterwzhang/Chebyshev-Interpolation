# initializeChebyshevApproximator ------------------------------------------------------------------
#' Initializes the Chebyshev Approximation
#'
#' @param D Dimensions of the Problem (integer)
#' @param N Highest Degree of the Polynomial (integer)
#' @param M Number of Interpolation Nodes in each dimension (integer)
#' @param bounds Bounds of the rectangle on which the function is approximated (list)
#' @param upper_b A vector of upper bounds (numeric)
#' @param lower_b A vector of lower bounds (numeric)
#' @return A list of the initialized approximation
#' @export

initializeChebyshevApproximator   <-    function(D, N, M = N + 1, bounds = NULL,
                                                upper_b = NULL, lower_b = NULL)
{
    stopifnot(M > N)
    if (all(is.null(upper_b), is.null(lower_b)))
    {
        upper_b   <-   unlist(lapply(bounds, function(ii) ii[2]))
        lower_b   <-   unlist(lapply(bounds, function(ii) ii[1]))
    }
    stopifnot(all(upper_b > lower_b))
    stopifnot(length(upper_b) == D & length(lower_b) == D)
    delta     <-   upper_b - lower_b

    # Calculate the Original M Chebyshev Nodes on [-1.1] and stack them
    nu   <-   unlist(lapply(1:M, function(i)
    {
        -1*cos((2*i-1)/(2*M) * base::pi)
    }))

    x_nodes   <-   lapply(1:D, function(d)
    {
        0.5*delta[d]*(nu + 1)+lower_b[d]
    })

    x_mat   <-   expand.grid(x_nodes, KEEP.OUT.ATTRS = TRUE)

    # Compute Tensor Product Basis and stack them
    T1   <-   calculateChebyshevPolynomials(nu, N)

    T2   <-   unlist(lapply(1:(N+1), function(i)
    {
        T1[,i] %*% T1[,i]
    }))

    B       <-   T1
    t_mat   <-   B
    if (D > 1)
    {
        for (d in 2:D)
        {
            t_mat   <-   kronecker(B, t_mat)  # Rows are stacked in T's rows
        }
    }

    a_mat   <-   solve(t(t_mat) %*% t_mat)%*% t(t_mat)
    theta   <-   matrix(0, nrow = (N+1)^D, ncol = 1)

    # Package the results into a list
    return(list(method   =   "Chebyshev-Approximation",
                D        =   D,        # Number of Dimensions
                N        =   N,        # Number of Degrees for Polynomial
                M        =   M,        # Number of Nodes
                bounds   =   bounds,   # Upper/Lower bounds in a list format
                upper_b  =   upper_b,  # Upper bound in a vector
                lower_b  =   lower_b,  # Lower bound in a vector
                delta    =   delta,    # Difference between upper and lower bound
                nu       =   nu,       # Vector of Approximation Nodes of length M
                X        =   x_mat,    # Stacked (row-wise) D-fold Cartesian Product matrix of nus of size M^D by D
                T        =   t_mat,    # Stacked (row-wise) Tensor-Basis Polynomials matrix of size M^D by (N+1)^D
                T1       =   T1,       # Pre-computed Chebyshev Polynomial Values of size M by (N+1)
                T2       =   T2,       # Pre-computed sum of squares values on nodes in a vector
                theta    =   theta,    # Chebyshev Coefficient Placeholder - all zeros
                A        =   a_mat     # Precomputed Weighting for calculating Coefficients of size (N+1)^D by M^D
                ))
}
