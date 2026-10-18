// Normal-gamma-gamma prior with heredity for pairwise interactions
// (Griffin and Brown, 2017, Example 5.1).
//
// beta_j ~ N(0, Psi_j / sd(X_j)^2)
// Psi_j = a1 * eta1_j for main effects
// Psi_jk = a2 * eta2_jk * h_jk for interactions, with
//   h_jk = eta1_j * eta1_k          (heredity = 1, strong)
//   h_jk = (eta1_j + eta1_k) / 2    (heredity = 2, weak)
// eta1_j ~ Gamma(lambda1, gamma1_j), gamma1_j ~ Gamma(c, (c - 1) / lambda1)
// eta2_jk ~ Gamma(lambda2, gamma2_jk), gamma2_jk ~ Gamma(c, (c - 1) / lambda2)
// a1 = lambda1 * d, a2 = lambda2 * d, lambda2 = r * lambda1
// lambda1 ~ Exponential(1), r ~ Beta(2, 6), p(d) proportional to (1 + d)^-2

data {
  int<lower=1> N;
  vector[N] y;
  int<lower=1> p;              // number of main effects
  int<lower=1> K;              // p + p * (p - 1) / 2
  matrix[N, K + 1] X;          // first column is the intercept
  int<lower=1> Ntest;
  vector[Ntest] ytest;
  matrix[Ntest, K + 1] Xtest;
  int<lower=1, upper=2> heredity;
  real<lower=1> c;             // tail parameter, c > 1 gives a finite mean
  real<lower=0> scale_sigma;
}
transformed data {
  int m = (p * (p - 1)) %/% 2;
  matrix[N, K] Xc;             // centred predictors without the intercept
  vector[K] means_X;
  vector[K] inv_sds_X;
  real ymean = mean(y);
  vector[N] yc = y - ymean;
  // Parents of each interaction, in the column order (1,2), (1,3), ..., (1,p), (2,3), ...
  array[m] int j_idx;
  array[m] int k_idx;
  {
    int q = 1;
    for (j in 1:(p - 1)) {
      for (k in (j + 1):p) {
        j_idx[q] = j;
        k_idx[q] = k;
        q += 1;
      }
    }
  }

  for (i in 1:K) {
    means_X[i] = mean(X[, i + 1]);
    inv_sds_X[i] = 1.0 / sd(X[, i + 1]);
    Xc[, i] = X[, i + 1] - means_X[i];
  }
}
parameters {
  real Intercept;              // intercept of the centred predictors
  vector[K] zbeta;
  real<lower=0> sigma;
  real<lower=0> lambda1;
  real<lower=0, upper=1> r;
  real<lower=0, upper=1> u_d;  // d = u_d / (1 - u_d)
  vector<lower=0>[p] eta1;
  vector<lower=0>[p] gamma1;
  vector<lower=0>[m] eta2;
  vector<lower=0>[m] gamma2;
}
transformed parameters {
  real d = u_d / (1.0 - u_d);
  real lambda2 = r * lambda1;
  real a1 = lambda1 * d;
  real a2 = lambda2 * d;
  vector[p] Psi_main;
  vector[m] Psi_int;
  vector[K] beta;

  for (j in 1:p) {
    Psi_main[j] = a1 * eta1[j];
  }
  for (q in 1:m) {
    real h;
    if (heredity == 1) {
      h = eta1[j_idx[q]] * eta1[k_idx[q]];
    } else {
      h = 0.5 * (eta1[j_idx[q]] + eta1[k_idx[q]]);
    }
    Psi_int[q] = a2 * eta2[q] * h;
  }
  for (j in 1:p) {
    beta[j] = sqrt(Psi_main[j]) * zbeta[j];
  }
  for (q in 1:m) {
    beta[p + q] = sqrt(Psi_int[q]) * zbeta[p + q];
  }
  beta = beta .* inv_sds_X;
}
model {
  target += normal_id_glm_lpdf(yc | Xc, Intercept, beta, sigma);
  target += std_normal_lpdf(zbeta);

  real d1 = (c - 1.0) / lambda1;
  target += gamma_lpdf(gamma1 | c, d1);
  target += gamma_lpdf(eta1 | lambda1, gamma1);

  real d2 = (c - 1.0) / (r * lambda1);
  target += gamma_lpdf(gamma2 | c, d2);
  target += gamma_lpdf(eta2 | (r * lambda1), gamma2);

  target += exponential_lpdf(lambda1 | 1.0);
  target += beta_lpdf(r | 2.0, 6.0);
  target += beta_lpdf(u_d | 1.0, 1.0);  // gives p(d) proportional to (1 + d)^-2
  target += student_t_lpdf(sigma | 3.0, 0.0, scale_sigma);
  target += normal_lpdf(Intercept | 0.0, 5.0);
}
generated quantities {
  real b_Intercept = ymean + Intercept - dot_product(means_X, beta);
  vector[N] mu = b_Intercept + X[, 2:(K + 1)] * beta;
  vector[Ntest] mu_test = b_Intercept + Xtest[, 2:(K + 1)] * beta;
  vector[N] log_lik;
  array[N] real y_tilde;
  vector[Ntest] log_lik_test;
  array[Ntest] real y_tilde_test;

  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | mu[n], sigma);
    y_tilde[n] = normal_rng(mu[n], sigma);
  }
  for (n in 1:Ntest) {
    log_lik_test[n] = normal_lpdf(ytest[n] | mu_test[n], sigma);
    y_tilde_test[n] = normal_rng(mu_test[n], sigma);
  }
}
