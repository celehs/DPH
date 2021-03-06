#' Accuracy Estimation
#' @param score risk score
#' @param X predictors
#' @param beta regression coefficients
#' @param lambda0 baseline hazards
#' @param w weights
#' @export
acc_est <- function(score, X, beta, lambda0, w = NULL) {
  eta <- c(as.matrix(X) %*% beta)  
  if (is.null(w)) w <- rep(1, length(eta))  
  K <- length(lambda0)
  Lambda0 <- cumsum(lambda0)  
  Surv0 <- exp(-Lambda0)  
  Surv <- outer(Surv0, exp(eta), "^")
  score_cut <- sort(unique(score), decreasing = TRUE)
  I <- 1 * outer(score, score_cut, ">=")
  FPR <- TPR <- matrix(NA, length(score_cut), K)
  AUC <- rep(NA, K)
  for (k in 1:K) {
    v1 <- w * Surv[k, ] / sum(w * Surv[k, ])
    v2 <- w * (1 - Surv[k, ]) / sum(w * (1 - Surv[k, ]))
    FPR[, k] <- colSums(v1 * I) 
    TPR[, k] <- colSums(v2 * I)
    AUC[k] <- sum(diff(c(0, FPR[, k])) * TPR[, k])
  }
  FPR <- tibble::as_tibble(data.frame(score_cut = score_cut, FPR))
  TPR <- tibble::as_tibble(data.frame(score_cut = score_cut, TPR))
  colnames(FPR)[-1] <- colnames(TPR)[-1] <- names(AUC) <- names(lambda0)
  list(FPR = FPR, TPR = TPR, AUC = AUC)
}
