
setwd("D:/Ohio State/SAMSI IMSM/EPA_SAMSI")

library(stats4)

likelihoodFunDensity <- function(phi, Y, dist, mu, sd){
  #Returns the likelihood function for the weighted average of densities
  #
  #Args:
  #  phi: Parameter to be estimated
  #  Y: AQS readings corresponding to DS estimates. Nx2 matrix.
  #  dist: Distance, measured in degrees of longitude, from each AQS point to
  #     boundaries of the overlap region. Nx2 matrix.
  #  mu: Downscaler estimates corresponding to AQS readings, Nx2 matrix.
  #  sd: Downscaler standard errors, Nx2 matrix.
  #
  #Returns:
  #  The log likelihood for phi
  w <- exp(-phi*dist)
  likelihood <- prod((1/(apply(w, 1, sum)))*(apply(w*dnorm(Y, mu, sd), 1, sum)))
  return(-(log(likelihood)))
}

likelihoodFunRV <-function(phi, Y, dist, mu, sd){  
  #Returns the likelihood function for the weighted average of densities
  #
  #Args:
  #  phi: Parameter to be estimated
  #  Y: AQS readings corresponding to DS estimates. Nx2 matrix.
  #  dist: Distance, measured in degrees of longitude, from each AQS point to
  #     boundaries of the overlap region. Nx2 matrix.
  #  mu: Downscaler estimates corresponding to AQS readings, Nx2 matrix.
  #  sd: Downscaler standard errors, Nx2 matrix.
  #
  #Returns:
  #  The log likelihood for phi
  w <- exp(-phi*dist)
  new.mu <- apply(w*mu, 1, sum)
  new.s <- apply((w^2)*(sd^2), 1, sum)
  new.s <- (new.s)^{1/2}
  likelihood <- prod(dnorm(Y, new.mu, new.s))
  return(-(log(likelihood)))
}

likelihoodFunRV2 <-function(a1, a2 ,dist.cen, Y, dist, mu, sd){  
  #Returns the likelihood function for the weighted average of densities
  #
  #Args:
  #  
  #  a1: Parameter to be estimated, intercept in phi
  #  a2: parameter to be estimated, slope in phi
  #  (phi: a1+a2xdist.cen)
  #  Y: AQS readings corresponding to DS estimates. Nx2 matrix.
  #  dist: Distance, measured in degrees of longitude, from each AQS point to
  #     boundaries of the overlap region. Nx2 matrix.
  #  mu: Downscaler estimates corresponding to AQS readings, Nx2 matrix.
  #  sd: Downscaler standard errors, Nx2 matrix.
  #  dist.cen:distnace measured from boundary center to the point. nx1 vector
  #Returns:
  #  The log likelihood for phi
  phi <- a1 + a2*dist.cen
  w <- exp(-phi*dist)
  new.mu <- apply(w*mu, 1, sum)
  new.s <- apply((w^2)*(sd^2), 1, sum)
  new.s <- (new.s)^{1/2}
  likelihood <- prod(dnorm(Y, new.mu, new.s))
  return(-(log(likelihood)))
}

smoothEstimate <- function(dist, mu, phi){
  #Returns combined estimates for each grid value in the intersection
  #
  #Args:
  #  dist: Distance, measured in degrees of longitude, from each grid point to
  #     boundaries of the overlap region. Nx2 matrix.
  #  mu: All Downscaler estimates, Nx2 matrix.
  #  phi: MLE of parameter phi
  #
  #Returns:
  # Combined estimates for each grid point, Nx1 vector
  
  w <- exp(-1*phi*dist)
  estimate <- (1/apply(w, 1, sum))*apply(w*mu, 1, sum)
  return(estimate)
}

smoothEstimate2 <- function(a1,a2, dist.cen, dist, mu){
  #Returns combined estimates for each grid value in the intersection
  #
  #Args:
  #  dist: Distance, measured in degrees of longitude, from each grid point to
  #     boundaries of the overlap region. Nx2 matrix.
  #  mu: All Downscaler estimates, Nx2 matrix.
  #  phi: MLE of parameter phi
  #
  #Returns:
  # Combined estimates for each grid point, Nx1 vector
  phi<-a1+a2*dist.cen
  w <- exp(-1*phi*dist)
  estimate <- (1/apply(w, 1, sum))*apply(w*mu, 1, sum)
  return(estimate)
}
