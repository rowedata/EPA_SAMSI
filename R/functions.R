
setwd("D:/Ohio State/SAMSI IMSM/EPA_SAMSI")

library(stats4)

likelihoodFunction <- function(phi, Y, d, mu, sd){
  #Returns the likelihood function for the weighted average of densities
  #
  #Args:
  #  phi: Parameter to be estimated
  #  Y: AQS readings corresponding to DS estimates. Nx2 matrix.
  #  d: Distance, measured in degrees of longitude, from each AQS point to
  #     boundaries of the overlap region. Nx2 matrix.
  #  mu: Downscaler estimates corresponding to AQS readings, Nx2 matrix.
  #  sd: Downscaler standard errors, Nx2 matrix.
  #
  #Returns:
  #  The log likelihood for phi
  w <- exp(-phi*d)
  R <- prod((1/(apply(w, 1, sum)))*(apply(w*dnorm(Y, mu, sd), 1, sum)))
  -(log(R))
}

smoothEstimate <- function(d, mu){
  #Returns combined estimates for each grid value in the intersection
  #
  #Args:
  #  d: Distance, measured in degrees of longitude, from each grid point to
  #     boundaries of the overlap region. Nx2 matrix.
  #  mu: All Downscaler estimates, Nx2 matrix.
  #
  #Returns:
  # Combined estimates for each grid point
  w <- exp(-1*phi.mle*d)
  f <- (1/apply(w, 1, sum))*apply(w*mu, 1, sum)
}
