// Regularized horseshoe prior for linear regression (Piironen and Vehtari, 2017).

functions {
  // Coefficients from the regularized horseshoe, see Appendix C.1 of
  // Piironen and Vehtari (2017), https://doi.org/10.1214/17-EJS1337SI
  vector horseshoe(vector z, vector lambda, real tau, real c2, vector inv_sds_X) {
    int K = rows(z);
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
    return z .* lambda_tilde * tau .* inv_sds_X;
  }
}
data {
  int<lower=1> N;
  vector[N] y;
  int<lower=1> K;              // number of coefficients
  matrix[N, K + 1] X;          // first column is the intercept
  int<lower=1> Ntest;
  vector[Ntest] ytest;
  matrix[Ntest, K + 1] Xtest;
  real<lower=0> hs_df;         // local degrees of freedom
  real<lower=0> hs_df_global;  // global degrees of freedom
  real<lower=0> hs_df_slab;    // slab degrees of freedom
  real<lower=0> hs_scale_slab; // slab scale
  real<lower=0> scale_sigma;
}
transformed data {
  real p0 = 0.5 * K;           // prior guess of the number of non-zero coefficients
  matrix[N, K] Xc;             // centred predictors without the intercept
  vector[K] means_X;
  vector[K] inv_sds_X;
  real ymean = mean(y);
  vector[N] yc = y - ymean;

  for (i in 1:K) {
    means_X[i] = mean(X[, i + 1]);
    inv_sds_X[i] = 1.0 / sd(X[, i + 1]);
    Xc[, i] = X[, i + 1] - means_X[i];
  }
}
parameters {
  real<lower=0> sigma;
  vector[K] zbeta;
  real Intercept;              // intercept of the centred predictors
  vector<lower=0>[K] hs_local;
  real<lower=0> hs_global;
  real<lower=0> hs_slab;
}
transformed parameters {
  vector[K] beta = horseshoe(zbeta, hs_local, hs_global, hs_scale_slab^2 * hs_slab, inv_sds_X);
}
model {
  target += normal_id_glm_lpdf(yc | Xc, Intercept, beta, sigma);
  target += std_normal_lpdf(zbeta);
  target += normal_lpdf(Intercept | 0, 5);
  hs_global ~ student_t(hs_df_global, 0, p0 / (K - p0) * sigma / sqrt(N));
  hs_slab ~ inv_gamma(0.5 * hs_df_slab, 0.5 * hs_df_slab);
  hs_local ~ student_t(hs_df, 0, 1);
  sigma ~ student_t(3, 0, scale_sigma);
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
