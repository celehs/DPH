acc_est <- function(score, lambda0) {
  K <- length(lambda0)
  Lambda0 <- cumsum(lambda0)  
  Surv0 <- exp(-Lambda0)  
  Surv <- outer(Surv0, exp(score), "^")
  score_cut <- sort(unique(score), decreasing = TRUE)
  I <- 1 * outer(score, score_cut, ">=")
  FPR <- TPR <- matrix(NA, length(score_cut), K)
  AUC <- rep(NA, K)
  for (k in 1:K) {
    w1 <- Surv[k, ] / sum(Surv[k, ])
    w2 <- (1 - Surv[k, ]) / sum(1 - Surv[k, ])
    FPR[, k] <- colSums(w1 * I) 
    TPR[, k] <- colSums(w2 * I)
    AUC[k] <- sum(diff(c(0, FPR[, k])) * TPR[, k])
  }
  return(AUC)
}
