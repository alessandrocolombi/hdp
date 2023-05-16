#' predictive_marginal
#'
#' @return [matrix] of size \code{n x length(grid)} containing the quantiles of level \code{0.025,0.5,0.975}.
#' @export
predictive <- function( idx_group, grid, fit,
                        priorMean, priorA, priorB, priorLambda, burnin = 0)
{

  n_iter <- length(fit$K) #number of iterations
  n_save <- n_iter - burnin #number of saved iterations to take into account. The first burnin must be discarded
  l_grid <- length(grid)  #length of the grid
  MIX    <- matrix(0, nrow=n_save, ncol=l_grid)

  # This loop computes the predictive distribution over a grid
  for(it in (burnin + 1):n_iter){

    # Get sampled values
    K_it    <- fit$K[it]                     # get the number of clusters
    mu_it   <- fit$mu[[it]]                  # get the mean, (mu_{1}^{(it)}, ..., mu_{M}^{(it)})
    sig2_it <- fit$sigma[[it]]^2               # get the variances, (sigma^2_{1}^{(it)}, ..., sigma^2_{M}^{(it)})
    q_it    <- fit$q_pred[[it]][idx_group,]  # get the weights of required group

    if(length(q_it)!=(K_it+1))
      stop("Error in predictive_marginal, length of q_it must be K_it + 1")

    # Kernel_grid is a  K_it x l_grid matrix, it contains the Normal kernels evauated over the grid
    # Kernel_grid[m,i] = Norm(grid[i] | mu_{m}^{(it)}, sigma^2_{m}^{(it)})
    Kernel_grid = t(sapply(1:K_it, simplify = "matrix",
                           function(m){
                             dnorm( x = grid, mean=mu_it[m], sd=sqrt(sig2_it[m]) )
                           }
    ))

    # Prior_grid is a vector of length l_grid, it contains the marginal prior evauated over the grid
    # Prior_grid[i] = nct(grid[i] | dof = nu0, loc = mu0, scale = sqrt(k0/(k0+1)*sigma0) ), where nct is the non-central student-t distribution

    sigma0 = priorB/priorA
    nu0 = 2*priorA
    scale = sqrt( (priorLambda+1)/(priorLambda) * sigma0 )
    if(scale <= 0)
      stop("The scale parameter has to be strictly positive.")

    Prior_grid = 1/scale * dt(x = (grid-priorMean)/scale, df = nu0 )

    # Compute predicted density at iteration it
    MIX[it-burnin,] <- q_it[K_it+1] * Prior_grid + q_it[1:K_it] %*% Kernel_grid
  }
  # Density estimation and credible bounds
  pred_est <- apply(MIX,2,quantile,prob=c(0.025,0.5,0.975))
  return(pred_est)

}


#' predictive_all_groups
#'
#' @return [list] of length \code{d} where each element is the return object of \code{\link{predictive}}.
#' @export
predictive_all_groups <- function( grid, fit,
                                   priorMean, priorA, priorB, priorLambda,
                                   burnin = 0)
{
  d = nrow(fit$q_pred[[1]])
  lapply(1:d, predictive, grid = grid, fit = fit,
                          priorMean = priorMean, priorA = priorA, priorB = priorB, priorLambda = priorLambda,
                          burnin = burnin)
}
