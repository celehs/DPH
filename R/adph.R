#' ADPH Log-likelihood
#' @param beta regression coefficients
#' @param lambda0 baseline hazards
#' @param sens sensitivities
#' @param spec specificities
#' @param t0 observed event time
#' @param d0 observed event indicator
#' @param X baseline covaraites
#' @export
adph_loglik <- function(beta, lambda0, sens, spec, t0, d0, X) {
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
  sum(log(lik))
}

#' Adjusted Discrete Proportional Hazards (ADPH) Model
#' @param time a vector of discrete time (e.g. 1, 2, 3, ...)
#' @param status a vector of event status (1 = observed, 0 = censored)
#' @param pred a vector/matrix of predictors (e.g. biomarkers)
#' @param sens sensitivity (a scalar or vector)
#' @param spec specificity (a scalar or vector)
#' @param sens_known indicator of whether sensitivity is known
#' @param spec_known indicator of whether specificity if known
#' @export
adph <- function(time, status, pred, sens, spec, sens_known = TRUE, spec_known = TRUE) {
  nt <- max(time)
  X <- as.matrix(pred)  
  np <- ncol(X)  
  if (length(sens) == 1) sens <- rep(sens, nt)
  if (length(spec) == 1) spec <- rep(spec, nt)
  fit0 <- dph(time, status, pred)
  start <- c(fit0$beta, stats::qlogis(fit0$lambda0))
  if (sens_known & spec_known) {
    # case 1: sens known, spec known
    fit <- stats::optim(
      par = start,
      fn = function(par) {
        (-1) * adph_loglik(
          beta = par[1:np],
          lambda0 = stats::plogis(par[np + 1:nt]),
          sens = sens,
          spec = spec, 
          t0 = c(time),
          d0 = c(status), 
          X = as.matrix(pred)
        )
      },
      method = "BFGS")      
  } else if (!sens_known & spec_known) {
    # case 2: sens unknown, spec known
    start <- c(start, stats::qlogis(sens))
    fit <- stats::optim(
      par = start,
      fn = function(par) {
        (-1) * adph_loglik(
          beta = par[1:np],
          lambda0 = stats::plogis(par[np + 1:nt]),
          sens = stats::plogis(par[np + nt + 1:nt]),
          spec = spec,
          t0 = c(time),
          d0 = c(status), 
          X = as.matrix(pred)
        )
      },
      method = "BFGS")     
    sens <- stats::plogis(fit$par[np + nt + 1:nt])  
  } else if (sens_known & !spec_known) {
    # case 3: sens known, spec unknown
    start <- c(start, stats::qlogis(spec))
    fit <- stats::optim(
      par = start,
      fn = function(par) {
        (-1) * adph_loglik(
          beta = par[1:np],
          lambda0 = stats::plogis(par[np + 1:nt]),
          sens = sens,
          spec = stats::plogis(par[np + nt + 1:nt]),
          t0 = c(time),
          d0 = c(status), 
          X = as.matrix(pred)
        )
      },
      method = "BFGS")     
    spec <- stats::plogis(fit$par[np + nt + 1:nt])  
  } else if (!sens_known & !spec_known) {
    # case 4: sens unknown, spec unknown    
    start <- c(start, stats::qlogis(sens), stats::qlogis(spec))
    fit <- stats::optim(
      par = start,
      fn = function(par) {
        (-1) * adph_loglik(
          beta = par[1:np],
          lambda0 = stats::plogis(par[np + 1:nt]),
          sens = stats::plogis(par[np + nt + 1:nt]),
          spec = stats::plogis(par[np + nt + nt + 1:nt]),
          t0 = c(time),
          d0 = c(status), 
          X = as.matrix(pred)
        )
      },
      method = "BFGS")     
    sens <- stats::plogis(fit$par[np + nt + 1:nt])     
    spec <- stats::plogis(fit$par[np + nt + nt + 1:nt])     
  }
  beta <- fit$par[1:np]
  lambda0 <- stats::plogis(fit$par[np + 1:nt])     
  acc <- acc_est(score = c(X %*% beta), X = X, beta = beta, lambda0 = lambda0)
  list(
    convergence = fit$convergence,
    beta = beta, 
    lambda0 = lambda0, 
    sens = sens, 
    spec = spec, 
    loglik = -fit$value, 
    FPR = acc$FPR,
    TPR = acc$TPR,
    AUC = acc$AUC
  )  
}
