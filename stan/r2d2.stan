// R2-D2 prior for linear regression (Zhang et al., 2022).
//
// R2 ~ Beta(R2_mean * R2_prec, (1 - R2_mean) * R2_prec)
// tau2 = sigma^2 * R2 / (1 - R2)
// phi ~ Dirichlet(R2_alpha)
// beta_j ~ N(0, phi_j * tau2 / sd(X_j)^2)

data {
  int<lower=1> N;
  vector[N] y;
  int<lower=1> K;              // number of coefficients
  matrix[N, K + 1] X;          // first column is the intercept
  int<lower=1> Ntest;
  vector[Ntest] ytest;
  matrix[Ntest, K + 1] Xtest;
  real<lower=0> R2_mean;
  real<lower=0> R2_prec;
  real<lower=0> R2_alpha;      // Dirichlet concentration
  real<lower=0> scale_sigma;
}
transformed data {
  matrix[N, K] Xc;             // centred predictors without the intercept
  vector[K] means_X;
  vector[K] inv_sds_X;
  vector[K] alpha = rep_vector(R2_alpha, K);
  real ymean = mean(y);
  vector[N] yc = y - ymean;

  for (i in 1:K) {
    means_X[i] = mean(X[, i + 1]);
    inv_sds_X[i] = 1.0 / sd(X[, i + 1]);
    Xc[, i] = X[, i + 1] - means_X[i];
  }
}
parameters {
  vector[K] zbeta;
  simplex[K] phi;
  real Intercept;              // intercept of the centred predictors
  real<lower=0, upper=1> R2;
  real<lower=0> sigma;
}
transformed parameters {
  real<lower=0> tau2 = sigma^2 * R2 / (1 - R2);
  vector[K] beta = zbeta .* sqrt(phi * tau2) .* inv_sds_X;
}
model {
  target += normal_id_glm_lpdf(yc | Xc, Intercept, beta, sigma);
  target += beta_lpdf(R2 | R2_mean * R2_prec, (1 - R2_mean) * R2_prec);
  target += dirichlet_lpdf(phi | alpha);
  target += std_normal_lpdf(zbeta);
  target += student_t_lpdf(sigma | 3, 0, scale_sigma);
  target += normal_lpdf(Intercept | 0, 5);
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
