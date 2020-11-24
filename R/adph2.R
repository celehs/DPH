#' ADPH2 Log-likelihood
#' @param beta regression coefficients
#' @param lambda0 baseline hazards
#' @param lambda0_s ...
#' @param t0 observed event time
#' @param d0 observed event indicator
#' @param X predictors
#' @param w weights
#' @export
adph2_loglik <- function(beta, lambda0, lambda0_s, t0, d0, X, w = NULL) {
  nt <- max(t0)
  n <- length(t0)  
  if (is.null(w)) w <- rep(1, n)
  eta <- c(as.matrix(X) %*% beta)
  O <- t(outer(1 - lambda0, exp(eta), "^"))
  S <- t(apply(O, 1, cumprod))  
  P <- cbind(1, S[, 1:(nt - 1)]) - S  
  A <- diag(lambda0_s[1], nt)
  for (k in 1:(nt - 1)) {
    for (kk in (k + 1):nt) {
      A[kk, k] <- prod((1 - lambda0_s)[1:(kk - k)]) * lambda0_s[kk - k + 1]
    }
  }
  P0 <- A %*% t(P)
  S0 <- 1 - apply(P0, 2, cumsum)
  P0 <- t(P0)
  S0 <- t(S0)  
  index <- cbind(1:n, t0)  
  loglik <- w * (d0 * log(P0[index]) + (1 - d0) * log(S0[index]))
  return(sum(loglik))
}

#' Adjusted Discrete Proportional Hazards (ADPH) Model
#' @param time a vector of discrete time (e.g. 1, 2, 3, ...)
#' @param status a vector of event status (1 = observed, 0 = censored)
#' @param pred a vector/matrix of predictors (e.g. biomarkers)
#' @param lambda0 baseline hazards
#' @param lambda0_s ...
#' @param n_ptb number of perturbations
#' @param seed random number generation seed
#' @export
adph2 <- function(time, status, pred, lambda0 = NULL, lambda0_s = NULL, 
                  n_ptb = 100, seed = 1) {
  n <- length(time)
  nt <- max(time)
  X <- as.matrix(pred)  
  np <- ncol(X)  
  fit0 <- dph(time, status, pred)
  if (is.null(lambda0)) lambda0 <- fit0$lambda0
  if (is.null(lambda0_s)) lambda0_s <- lambda0
  start <- c(fit0$beta, stats::qlogis(lambda0), stats::qlogis(lambda0_s)) 
  fit <- stats::optim(
    par = start,
    fn = function(par) {
      (-1) * adph2_loglik(
        beta = par[1:np],
        lambda0 = stats::plogis(par[np + 1:nt]),
        lambda0_s = stats::plogis(par[np + nt + 1:nt]),
        t0 = c(time),
        d0 = c(status), 
        X = X
      )
    },
    method = "BFGS"
  )     
  beta <- fit$par[1:np]
  lambda0 <- stats::plogis(fit$par[np + 1:nt]) 
  lambda0_s <- stats::plogis(fit$par[np + nt + 1:nt])
  loglik <- -fit$value
  AUC <- acc_est(score = c(X %*% beta), lambda0 = lambda0)
  names(AUC) <- names(lambda0) 
  beta_ptb <- data.frame(matrix(NA, n_ptb, np))
  lambda0_ptb <- data.frame(matrix(NA, n_ptb, nt))
  lambda0_s_ptb <- data.frame(matrix(NA, n_ptb, nt))
  colnames(beta_ptb) <- names(beta)
  colnames(lambda0_ptb) <- names(lambda0)
  colnames(lambda0_s_ptb) <- names(lambda0_s)
  set.seed(seed)
  for (i in 1:n_ptb) {
    w <- rexp(n)
    start <- c(beta, stats::qlogis(lambda0), stats::qlogis(lambda0_s)) 
    fit <- stats::optim(
      par = start,
      fn = function(par) {
        (-1) * adph2_loglik(
          beta = par[1:np],
          lambda0 = stats::plogis(par[np + 1:nt]),
          lambda0_s = stats::plogis(par[np + nt + 1:nt]),
          t0 = c(time),
          d0 = c(status), 
          X = X,
          w = w
        )
      },
      method = "BFGS"
    )
    beta_ptb[i, ] <- fit$par[1:np]
    lambda0_ptb[i, ] <- stats::plogis(fit$par[np + 1:nt]) 
    lambda0_s_ptb[i, ] <- stats::plogis(fit$par[np + nt + 1:nt])    
  }
  return(list(
    beta = beta, 
    lambda0 = lambda0, 
    lambda0_s = lambda0_s, 
    loglik = loglik, 
    AUC = AUC,
    beta_ptb = tibble::as_tibble(beta_ptb),
    lambda0_ptb = tibble::as_tibble(lambda0_ptb),
    lambda0_s_ptb = tibble::as_tibble(lambda0_s_ptb)
  ))  
}
