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

adph2_fit <- function(time, status, X, w, start, method = "BFGS") {
  nt <- max(time)  
  np <- ncol(X) 
  obj <- stats::optim(
    par = start,
    fn = function(par) {
      - adph2_loglik(
        beta = par[1:np],
        lambda0 = stats::plogis(par[np + 1:nt]),
        lambda0_s = stats::plogis(par[np + nt + 1:nt]),
        t0 = c(time),
        d0 = c(status), 
        X = X,
        w = w
      )
    },
    method = method
  )     
  loglik <- -obj$value
  beta <- obj$par[1:np]
  lambda0 <- stats::plogis(obj$par[np + 1:nt]) 
  lambda0_s <- stats::plogis(obj$par[np + nt + 1:nt])
  AUC <- auc(score = c(X %*% beta), X = X, beta = beta, lambda0 = lambda0)
  names(AUC) <- names(lambda0)   
  return(list(
    loglik = loglik,
    beta = beta, 
    lambda0 = lambda0,
    lambda0_s = lambda0_s,
    AUC = AUC))
}

#' Adjusted Discrete Proportional Hazards (ADPH) Model
#' @param time a vector of discrete time (e.g. 1, 2, 3, ...)
#' @param status a vector of event status (1 = observed, 0 = censored)
#' @param pred a vector/matrix of predictors (e.g. biomarkers)
#' @param lambda0 baseline hazards
#' @param lambda0_s ...
#' @param n_ptb number of perturbations
#' @param seed random number generation seed
#' @param show_progress TRUE or FALSE (default)
#' @export
adph2 <- function(time, status, pred, lambda0 = NULL, lambda0_s = NULL, 
                  n_ptb = 200, seed = 1, show_progress = FALSE) {
  n <- length(time)
  nt <- max(time)
  X <- as.matrix(pred)  
  np <- ncol(X)  
  fit0 <- dph(time, status, pred)
  if (is.null(lambda0)) lambda0 <- fit0$lambda0
  if (is.null(lambda0_s)) lambda0_s <- lambda0
  start <- c(fit0$beta, stats::qlogis(lambda0), stats::qlogis(lambda0_s)) 
  fit <- adph2_fit(
    time = time,
    status = status,
    X = X,
    w = rep(1, n),
    start = start
  )
  ptb <- NULL
  ptb$beta <- data.frame(matrix(NA, n_ptb, np))
  ptb$lambda0 <- data.frame(matrix(NA, n_ptb, nt))
  ptb$lambda0_s <- data.frame(matrix(NA, n_ptb, nt))
  ptb$AUC <- data.frame(matrix(NA, n_ptb, nt))
  colnames(ptb$beta) <- names(fit$beta)
  colnames(ptb$lambda0) <- names(fit$lambda0)
  colnames(ptb$lambda0_s) <- names(fit$lambda0_s)
  colnames(ptb$AUC) <- names(fit$AUC)
  set.seed(seed)
  start <- c(fit$beta, fit$lambda0, fit$lambda0_s) 
  for (i in 1:n_ptb) {
    if (show_progress) print(paste0("perturbation: ", i, "/", n_ptb))
    fit_ptb <- adph2_fit(
      time = time,
      status = status,
      X = X,
      w = stats::rexp(n),
      start = start
    ) 
    ptb$beta[i, ] <- fit_ptb$beta
    ptb$lambda0[i, ] <- fit_ptb$lambda0
    ptb$lambda0_s[i, ] <- fit_ptb$lambda0_s   
    ptb$AUC[i, ] <- fit_ptb$AUC
  }
  ptb$beta <- tibble::as_tibble(ptb$beta)
  ptb$lambda0 <- tibble::as_tibble(ptb$lambda0)
  ptb$lambda0_s <- tibble::as_tibble(ptb$lambda0_s)
  ptb$AUC <- tibble::as_tibble(ptb$AUC)
  return(list(fit = fit, ptb = ptb))  
}
