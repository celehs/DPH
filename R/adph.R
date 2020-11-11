#' ADPH Log-likelihood
#' @param lambda0 baseline hazards
#' @param beta regression coefficients
#' @param sens sensitivities
#' @param spec specificities
#' @param t0 observed event time
#' @param d0 observed event indicator
#' @param X baseline covaraites
#' @export
adph_loglik <- function(lambda0, beta, sens, spec, t0, d0, X) {
  theta <- sens
  phi <- spec
  n <- length(t0)
  eta <- c(as.matrix(X) %*% beta)
  gamma0 <- stats::qlogis(lambda0)
  M <- t(outer(1 + exp(gamma0), -exp(eta), "^"))
  lik <- rep(NA, n)
  for (i in 1:n) {
    phi0 <- phi[t0[i]]
    Gamma <- prod(phi[1:t0[i]]) *
      ((1 - phi0) / phi0)^d0[i]
    theta0 <- theta[t0[i]]
    Delta <- prod(1 - theta[1:t0[i]]) *
      (theta0 / (1 - theta0))^d0[i] # k = 1
    lik[i] <- prod(M[i, 1:t0[i]]) * Gamma
    lik[i] <- lik[i] + (1 - M[i, 1]) * Delta
    if (t0[i] > 1) {
      for (k in 2:t0[i]) {
        theta0 <- theta[t0[i] - k + 1]
        Delta <- prod(phi[1:(k - 1)]) *
          prod(1 - theta[0:(t0[i] - k) + 1]) *
          (theta0 / (1 - theta0))^d0[i]
        lik[i] <- lik[i] +
          prod(M[i, 1:(k - 1)]) *
          (1 - M[i, k]) * Delta
      }
    }
  }
  return(sum(log(lik)))
}

#' Adjusted Discrete Proportional Hazards (ADPH) Model
#' @param time a vector of discrete time (e.g. 1, 2, 3, ...)
#' @param status a vector of event status (1 = observed, 0 = censored)
#' @param pred a vector/matrix of predictors (e.g. biomarkers)
#' @param sens sensitivity (a scalar or vector)
#' @param spec specificity (a scalar or vector)
#' @export
adph <- function(time, status, pred, sens, spec) {
  t_max <- max(time)
  if (length(sens) == 1) sens <- rep(sens, t_max)
  if (length(spec) == 1) spec <- rep(spec, t_max)
  fit0 <- dph(time, status, pred)
  start <- c(stats::qlogis(fit0$lambda0), fit0$beta)
  fit <- stats::optim(
    par = start,
    fn = function(par) {
      - adph_loglik(
        lambda0 = stats::plogis(par[1:t_max]),
        beta = par[-(1:t_max)],
        sens = sens,
        spec = spec, 
        t0 = c(time),
        d0 = c(status), 
        X = as.matrix(pred)
      )
    },
    method = "BFGS")  
  lambda0 <- stats::plogis(fit$par[1:t_max])
  beta <- fit$par[-(1:t_max)]
  list(lambda0 = lambda0, beta = beta, sens = sens, spec = spec)  
}
