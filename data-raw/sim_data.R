mkdata <- function(n, beta, lambda0, lambda0_s, seed = 1) {
  set.seed(seed)
  z <- rnorm(n, 2, 1) 
  X <- matrix(z)
  eta <- c(X %*% beta)
  O <- t(outer(1 - lambda0, exp(eta), "^"))
  S <- t(apply(O, 1, cumprod))
  P <- cbind(1, S[, 1:4]) - S
  # range(apply(P, 1, sum))
  A <- diag(lambda0_s[1], 5)
  for (k in 1:4) {
    for (kk in (k + 1):5) {
      A[kk, k] <- prod((1 - lambda0_s)[1:(kk - k)]) * lambda0_s[kk - k + 1]
    }
  }
  time0 <- rep(NA, n)
  for (k in 1:n) {
    p <- P[k, ]
    p0 <- A %*% matrix(p, ncol = 1)
    M <- expand.grid(l = 1:6, w = NA)
    M[, "w"] <- c(p0, 1 - sum(p0))
    idx <- sample(1:nrow(M), size = 1, prob = M[, "w"])
    time0[k] <- M[idx, "l"]
  }
  c <- 5 # censoring time
  d0 <- 1 * (time0 <= c)
  t0 <- apply(cbind(time0, c), 1, min)
  # print(addmargins(table(d0, t0)))
  DF <- data.frame(
    time = as.integer(t0), 
    status = as.integer(d0), 
    pred1 = z, 
    pred2 = rbinom(n, 1, 0.5))
  tibble::as_tibble(DF)
}

sim_data <- mkdata(
  n = 1000, 
  beta = 1, 
  lambda0 = 0.05 * (1:5), 
  lambda0_s = 0.1 * (8:4))

usethis::use_data(sim_data, overwrite = TRUE)
