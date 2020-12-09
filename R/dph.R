#' Standard Discrete Proportional Hazards Model
#' @param time a vector of discrete time (e.g. 1, 2, 3, ...)
#' @param status a vector of event status (1 = observed, 0 = censored)
#' @param pred a vector/matrix of predictors (e.g. biomarkers)
#' @export
dph <- function(time, status, pred) {
  t <- c(time)
  d <- c(status)
  X <- as.matrix(pred)
  n <- nrow(X)
  t_max <- max(t)
  stopifnot(t > 0)
  # stopifnot(d %in% c(0, 1))
  stopifnot(length(t) == n)
  stopifnot(length(d) == n)  
  L <- vector("list", n)
  # create pseudo data
  for (i in 1:n) {
    y <- rep(0, t[i])
    if (d[i] == 1) y[t[i]] <- 1
    L[[i]] <- cbind(y, z = 1:t[i])
  }
  DF <- cbind(do.call("rbind", L), X[rep(seq(t), t), ])
  y <- DF[, 1]
  z <- factor(DF[, 2])
  Z <- DF[, -(1:2)]
  fit <- stats::glm(y ~ z + Z - 1, family = stats::binomial("cloglog"))
  alpha <- fit$coef[1:t_max]
  beta <- fit$coef[-(1:t_max)]
  lambda0 <- 1 - exp(-exp(alpha))
  names(lambda0) <- paste0("time", 1:t_max)
  names(beta) <- colnames(Z)
  acc <- acc_est(score = c(X %*% beta), X = X, beta = beta, lambda0 = lambda0)
  return(list(
    beta = beta, 
    lambda0 = lambda0, 
    loglik = c(stats::logLik(fit)),
    FPR = acc$FPR,
    TPR = acc$TPR,
    AUC = acc$AUC
  ))
}
