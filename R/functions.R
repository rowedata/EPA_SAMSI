
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
  #Returns the likelihood function for the weighted average of random variables
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
  w <- (1/apply(w, 1, sum))*w
  new.mu <- apply(w*mu, 1, sum)
  new.s <- apply((w^2)*(sd^2), 1, sum)
  new.s <- (new.s)^{1/2}
  likelihood <- prod(dnorm(Y, new.mu, new.s))
  return(-(log(likelihood)))
}

likelihoodFunRV2 <-function(a1, a2 ,dist.cen, Y, dist, mu, sd){  
  #Returns the likelihood function for the weighted average of random variables
  #  accounting for edges of interesection
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
  #  dist.cen: distnace measured from intersection center to the point. nx1 vector
  #
  #Returns:
  #  The log likelihood for phi
  phi <- a1 + a2*dist.cen
  w <- exp(-phi*dist)
  w <- (1/apply(w, 1, sum))*w
  new.mu <- apply(w*mu, 1, sum)
  new.s <- apply((w^2)*(sd^2), 1, sum)
  new.s <- (new.s)^{1/2}
  likelihood <- prod(dnorm(Y, new.mu, new.s))
  return(-(log(likelihood)))
}

likelihoodFunRV3 <-function(beta.0, beta.1 , alpha.0, alpha.1, d.1.centre, 
                            d.2.centre, d.1, d.2, Y, mu, sd){  
  #Returns the likelihood function for the weighted average of random variables
  #  accounting for edges of interesection
  #
  #Args:
  #  
  #  beta.0: Intercept term for phi
  #  beta.1: Slope for d.1.centre
  #  d.1.centre: Latitudinal distance to centre of intersection
  #  alpha.0: Intercept term for theta
  #  alpha.1: Slope for d.2.centre
  #  d.2.centre: Longitudinal distance to centre of intersection
  #  d.1: Latitudinal distances to centre of region 1 and 2, Nx2 matrix
  #  d.2: Longitudinal distances to centre of region 1 and 2, Nx2 matrix
  #  Y: AQS readings corresponding to DS estimates
  #  mu: Downscaler estimates corresponding to AQS readings, Nx2 matrix
  #  sd: Downscaler standard errors, Nx2 matrix
  #Returns
  #The log likelihood for parameters
  phi <- beta.0 + beta.1*d.1.centre
  theta <- alpha.0 + alpha.1*d.2.centre
  w <- exp(-(phi*d.1 + theta*d.2))
  w <- (1/apply(w, 1, sum))*w
  new.mu <- apply(w*mu, 1, sum)
  new.s <- sqrt(apply((w^2)*(sd^2), 1, sum))
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

smoothEstimate2 <- function(a1, a2, dist.cen, dist, mu){
  #Returns combined estimates for each grid value in the intersection
  #
  #Args:
  #  a1: Intercept for phi
  #  a2: Slope for phi vs. distance from centre of intersection
  #  dist: Distance, measured in degrees of longitude, from each grid point to
  #     boundaries of the overlap region. Nx2 matrix.
  #  mu: All Downscaler estimates, Nx2 matrix.
  #  phi: MLE of parameter phi
  #  dist.cen: Distance from centre of intersection to DS estimate
  #
  #Returns:
  # Combined estimates for each grid point, Nx1 vector
  phi <- a1 + a2*dist.cen
  w <- exp(-1*phi*dist)
  estimate <- (1/apply(w, 1, sum))*apply(w*mu, 1, sum)
  return(estimate)
}

smoothEstimate3 <- function(beta.0, beta.1 , alpha.0, alpha.1, d.1.centre, 
                            d.2.centre, d.1, d.2, mu){
  #Returns the likelihood function for the weighted average of random variables
  #  accounting for edges of interesection
  #
  #Args:
  #  
  #  beta.0: Intercept term for phi
  #  beta.1: Slope for d.1.centre
  #  d.1.centre: Latitudinal distance to centre of intersection
  #  alpha.0: Intercept term for theta
  #  alpha.1: Slope for d.2.centre
  #  d.2.centre: Longitudinal distance to centre of intersection
  #  d.1: Latitudinal distances to centre of region 1 and 2, Nx2 matrix
  #  d.2: Longitudinal distances to centre of region 1 and 2, Nx2 matrix
  #  mu: Downscaler estimates corresponding to AQS readings, Nx2 matrix
  #Returns
  #Spliced DS estimates
  phi <- beta.0 + beta.1*d.1.centre
  theta <- alpha.0 + alpha.1*d.2.centre
  w <- exp(-(phi*d.1 + theta*d.2))
  estimate <- (1/apply(w, 1, sum))*apply(w*mu, 1, sum)
}

standardErrorsModel3 <- function(sd, beta.0, beta.1, dist.cen, dist){
  #Returns standard errors for estimates from model 3
  #
  #Args:
  #  a1: Intercept for phi
  #  a2: Slope for phi vs. distance from centre of intersection
  #  dist: Distance, measured in degrees of longitude, from each grid point to
  #     boundaries of the overlap region. Nx2 matrix.
  #  mu: All Downscaler estimates, Nx2 matrix.
  #  phi: MLE of parameter phi
  #  dist.cen: Distance from centre of intersection to DS estimate
  #
  #Returns:
  # Combined estimates for each grid point, Nx1 vector
  phi <- beta.0 + beta.1*dist.cen
  w <- exp(-phi*dist)
  se <- sqrt((1/apply((w^2), 1, sum))*apply((w^2)*(sd^2), 1, sum))
  return(se)
}
