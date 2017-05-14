#' Create random sample from {-1, 1}
#'
#' Uses stats::rbinom to create a random {-1, 1} vector.
#'
#' @param n The number of samples.
#' @return  An n-length vector from {-1, 1}.
#' @export
rsign <- function(n){
  bin <- rbinom(n, 1, 0.5)  
  2 * bin - rep(1, n)
}

#' Jitter the parameters of a time series model
#'
#' @param mod A time series model.
#' @param percent Max percentage change in parameter values. 
#' @return A times series with randonly perturbed. 
#'         parameters in mod 
#' @export
jitter_params <- function(mod, percent) UseMethod("jitter_params")


#' @export
jitter_params.mackeyglass <- function(mod, percent = 0.10){
  len <- length(mod)
  # subtract or add (based on rsign) the perturbed parameter
  delta <- percent*runif(length(mod),0,1) * rsign(length(mod))
  vmod <- unlist(mod)*delta + unlist(mod)
  names(vmod) <- NULL 
  mackeyglass(tau = floor(vmod[1]), 
              beta = vmod[2], 
              gamma = vmod[3], 
              n = vmod[4], 
              noise = vmod[5], 
              init =vmod[6])
}


#' @export 
jitter_params.farima <- function(mod, percent = 0.10){
  ar <- mod$ar + percent*runif(length(mod$ar), 0, 1) * mod$ar * rsign(length(mod$ar))
  ma <- mod$ma + percent*runif(length(mod$ma), 0, 1) * mod$ma * rsign(length(mod$ma))
  d  <- mod$d  + percent*runif(1, 0, 1)*mod$d
  farima(ar, ma, d) 
}

#' @export 
jitter_params.arma <- function(mod, percent = 0.10){
  ar <- mod$ar + percent*runif(length(mod$ar), 0, 1) * mod$ar * rsign(length(mod$ar))
  ma <- mod$ma + percent*runif(length(mod$ma), 0, 1) * mod$ma * rsign(length(mod$ma))
  arma(ar, ma) 
}


#' @export 
jitter_params.logistic <- function(mod, percent = 0.01){
  delta <- percent*runif(1, 0, 1)*unlist(mod)
  r <- mod$r + delta
  if(r >= 4)  r <- mod$r - delta
  logistic(r) 
}

#' @export 
jitter_params.weierstrass <- function(mod, percent = 0.10){
  delta <- percent*runif(1, 0, 1) * mod$a * rsign(1)
  a <- mod$a + delta
  weierstrass(a) 
}

#' @export 
jitter_params.default <- function(mod, percent = 0.10){
  mod
}
