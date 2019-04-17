
# evaluateChebyshev -------------------------------------------------------------------------------
#' Evaluates the Chebyshev Approximation for a matrix (or a vector) of points
#'
#' Option for parallelized evaluation for many points to evaluate
#'
#' @param x Points to evaluate with size Points by Dimensions (matrix)
#' @param cheb List of item from initalizeChebyshevApproximator (list)
#' @param parallel Boolean flag for parallelization (logical)
#' @param numcores Cores for parallelization (integer)
#' @return A vector of predictions for each point of x
#' @export

evaluateChebyshev   <-   function(x, cheb, parallel = FALSE, numcores = 1L)
{
    stopifnot(class(x) == "matrix")
    stopifnot(ncol(x)  == cheb$D)
    xi   <-   sapply(1:cheb$D, function(d)
    {
        2*((x[,d]-cheb$lower_b[d])/cheb$delta[d]) - 1.0
    })

    if(class(xi) != "matrix")
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

    y   <-   NULL
    if (parallel)
    {
        y   <-   unlist(mclapply(1:nrow(x), function(k)
        {
            t1   <-   Ti[[1]][k,]
            if (cheb$D >= 2)
            {
                for (d in 2:cheb$D)
                {
                    t1   <-   kronecker(Ti[[d]][k,], t1)
                }
            }
            t1 %*% cheb$theta
        }, mc.cores = numcores))
    } else
    {
        y   <-   unlist(lapply(1:nrow(x), function(k)
        {
            t1   <-   Ti[[1]][k,]
            if (cheb$D >= 2)
            {
                for (d in 2:cheb$D)
                {
                    t1   <-   kronecker(Ti[[d]][k,], t1)
                }
            }
            t1 %*% cheb$theta
        }))
    }

    return(y)
}

# -------------------------------------------------------------------------------------------------
