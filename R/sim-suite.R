#' Generate A Default suite of function models.
#'
#' @return A list of two lists of function models with
#'          the same functions and varying default parameters.
#' @export
sim_suite <- function(){
  # The time serires models : arma, arma, farima, logistic, quadratic, Mackey-Glass

  farima1 <- farima(ar = c(0.1, -0.5), ma = c(0.6, 0.01), d = 0.3) 
  farima2 <- farima(ar = c(0.2, -0.4), ma = c(0.4, 0.02), d = 0.3)  

  arma1 <- arma(ar = c(-0.1, 0.3), ma = c(0.2, 0.1))
  arma2 <- arma(ar = c(0.4, 0.3), ma = c(0.1, -0.5))

  arma3 <- arma(ar  = c(0.5, -0.6, 0.9), ma = c(0.5, 0.6))
  arma4 <- arma(ar = c(-0.2, -0.3, -0.8), ma = c(0.2, 0.1))

  log1 <- logistic(r = 3.87)
  log2 <- logistic(r = 3.70)

  cauchy1 <- cauchy(alpha = 0.3, beta = 0.5) 
  cauchy2 <- cauchy(alpha = 0.4, beta = 0.7)
  # mg1 <- mackeyglass(tau = 27, 
  #                    beta = 0.3, 
  #                    gamma = 0.1,
  #                    n = 9, 
  #                    init = rep(0.2),
  #                    noise = 1)

  # mg2 <- mackeyglass(tau = 17, 
  #                    beta = 0.35, 
  #                    gamma = 0.2,
  #                    n = 8, 
  #                    init = rep(0.2),
  #                    noise = 1)

  weier1 <- weierstrass(a = 0.4, b = 4)
  weier2 <- weierstrass(a = 0.8, b = 4)

  fbm1 <- fBm(H = 0.1)
  fbm2 <- fBm(H = 0.3) 

  group1 <- list(ARMA1 = arma1, 
                 # arma2 = arma2, k
                 Logistic1     = log1, 
                 Weierstrass1  = weier1,
                 cauchy1       = cauchy1, 
                 # Mackey_Glass1 = mg1, 
                 FARIMA1       = farima1, 
                 fBm1          = fbm1)
  group2 <- list(ARMA2 = arma2, 
                 # arma2 = arma2, k
                 Logistic2     = log2, 
                 Weierstrass2  = weier2,
                 Cauchy1       = cauchy1,
                 # Mackey_Glass2 = mg2, 
                 FARIMA2       = farima2,
                 fBm2          = fbm2)



  return(list(group1 = group1, group2 = group2))
}
#' Returns the names corresonding to functions of \code{sim_suite()}
#'
#'
#' @return  A vector of function names
#' @export
sim_names <- function() {
  c("ARMA", "Logistic", "Weierstrass", "Cauchy", "FARIMA", "fBm")
}