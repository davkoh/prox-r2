# Data-generating processes of the simulation study.

# Stationary AR(1) covariance with unit innovation variance.
ar1_cov <- function(p, rho) toeplitz(rho^(0:(p - 1))) / (1 - rho^2)

gen_X <- function(n, p, rho, seed) {
  set.seed(seed)
  X <- mvtnorm::rmvnorm(n, sigma = ar1_cov(p, rho))
  colnames(X) <- paste0("x", 1:p)
  X
}

# Pairwise products in the column order (1,2), (1,3), ..., (1,p), (2,3), ...
gen_X_int <- function(X) {
  pairs <- combn(ncol(X), 2)
  X_int <- X[, pairs[1, ], drop = FALSE] * X[, pairs[2, ], drop = FALSE]
  colnames(X_int) <- paste(colnames(X)[pairs[1, ]], colnames(X)[pairs[2, ]], sep = "_")
  X_int
}

# Main effects are dense (all 0.05) or sparse (0.5 for the first and last
# three, 0 otherwise). Interactions satisfy weak heredity ("correct") or
# violate it ("misp").
gen_beta <- function(p, dgp_type, seed) {
  set.seed(seed)
  b_main <- if (startsWith(dgp_type, "dense")) {
    rep(0.05, p)
  } else {
    c(rep(0.5, 3), rep(0, p - 6), rep(0.5, 3))
  }
  names(b_main) <- paste0("x", 1:p)

  pairs <- combn(p, 2)
  b_j <- b_main[pairs[1, ]]
  b_k <- b_main[pairs[2, ]]
  bound <- 0.5 * (abs(b_j) + abs(b_k))

  b_int <- if (endsWith(dgp_type, "correct")) {
    # |b_jk| <= (|b_j| + |b_k|) / 2
    out <- numeric(length(bound))
    for (i in which(bound != 0)) {
      out[i] <- truncnorm::rtruncnorm(1, a = -bound[i], b = bound[i], mean = 0, sd = 3)
    }
    out
  } else {
    # b_jk equals the larger parent if one parent is zero, else the bound
    unname(ifelse(b_j == 0 | b_k == 0, pmax(abs(b_j), abs(b_k)), bound))
  }
  names(b_int) <- paste0("x", pairs[1, ], "_x", pairs[2, ])

  list(main = b_main, int = b_int)
}

# One replicate with n training and n test observations. sigma gives the
# target R2 = b'X'Xb / (b'X'Xb + n sigma^2) on the training data.
simulate_replicate <- function(n, p, rho, r2, dgp_type, seed) {
  X <- gen_X(n, p, rho, seed)
  X_int <- gen_X_int(X)
  beta <- gen_beta(p, dgp_type, seed)
  X_test <- gen_X(n, p, rho, seed + 1)
  X_int_test <- gen_X_int(X_test)

  X_full <- cbind(X, X_int)
  b <- c(beta$main, beta$int)
  sigma <- sqrt(drop(t(X_full %*% b) %*% X_full %*% b) / n * (1 - r2) / r2)

  # The noise draws restart the random number stream at seed + 1.
  set.seed(seed + 1)
  y <- as.numeric(X %*% beta$main + X_int %*% beta$int + rnorm(n, 0, sigma))
  y_test <- as.numeric(X_test %*% beta$main + X_int_test %*% beta$int + rnorm(n, 0, sigma))

  list(X = X, X_int = X_int, y = y,
       X_test = X_test, X_int_test = X_int_test, y_test = y_test,
       beta_main = beta$main, beta_int = beta$int, sigma = sigma)
}
