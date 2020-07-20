# evaluateChebyshev_T -----------------------------------------------------------------------------
#' Evaluates the Chebyshev Approximation for a matrix (or a vector) of points and returns
#' the underlying basis function values instead of the interpolation values
#'
#' Option for parallelized evaluation for many points to evaluate
#'
#' @param x Points to evaluate with size Points by Dimensions (matrix)
#' @param cheb List of item from initalizeChebyshevApproximator (list)
#' @param parallel Boolean flag for parallelization (logical)
#' @param numcores Cores for parallelization (integer)
#' @return A matrix of the underlying basis function values
#' @export

evaluateChebyshev_T   <-   function(x, cheb, parallel = FALSE, numcores = 1L)
{
    stopifnot(any(class(x) == "matrix"))
    stopifnot(ncol(x)  == cheb$D)

    xi   <-   sapply(1:cheb$D, function(d)
    {
        2*((x[,d]-cheb$lower_b[d])/cheb$delta[d]) - 1.0
    })

    if(!any(class(xi) == "matrix"))
    {
        if (nrow(xi) == 1)
        {
            # One Point Case
            xi   <-   matrix(xi, ncol = cheb$D)
        } else
        {
            # One Dimension Case
            xi   <-   matrix(xi, nrow = cheb$D)
        }
    }


    Ti   <-   lapply(1:cheb$D, function(d)
    {
        calculateChebyshevPolynomials(xi[,d], cheb$N)
    })

    T   <-   NULL
    if (parallel)
    {
        T   <-   do.call(rbind, mclapply(1:nrow(x), function(k)
        {
            t1   <-   Ti[[1]][k,]
            if (cheb$D >= 2)
            {
                for (d in 2:cheb$D)
                {
                    t1   <-   kronecker(Ti[[d]][k,], t1)
                }
            }
            t1
        }, mc.cores = numcores))
    } else
    {
        T   <-   do.call(rbind, lapply(1:nrow(x), function(k)
        {
            t1   <-   Ti[[1]][k,]
            if (cheb$D >= 2)
            {
                for (d in 2:cheb$D)
                {
                    t1   <-   kronecker(Ti[[d]][k,], t1)
                }
            }
            t1
        }))
    }


    return(T)
}

# -------------------------------------------------------------------------------------------------
