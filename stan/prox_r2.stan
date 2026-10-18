// Prox-R2 prior for linear regression with pairwise interactions.
//
// R2 ~ Beta(R2_mean * R2_prec, (1 - R2_mean) * R2_prec)
// tau2 = sigma^2 * R2 / (1 - R2)
// z ~ N(0, I) with the penalty -||z - Pi_C(z)||^2 / (2 * proximal_lambda)
// phi = softmax(Pi_C(z))
// beta_j ~ N(0, phi_j * tau2 / s_j^2)
//
// Pi_C is the projection onto the weak heredity set
// C = {x : 2 w_jk x_jk <= w_j x_j + w_k x_k for all j < k}, with w = 1 / s^2.
// Dykstra's algorithm computes the projection.

functions {
  // Position of interaction (j, k), j < k, in the coefficient vector.
  int pair_index(int j, int k, int p) {
    return p + ((k - 1) * (k - 2)) %/% 2 + j;
  }

  vector dykstra(vector g, vector w, int p, int max_iter) {
    int m = (p * (p - 1)) %/% 2;
    vector[num_elements(g)] x = g;
    // Increments of the previous cycle (three non-zero entries per constraint)
    array[m] real Ij = rep_array(0.0, m);
    array[m] real Ik = rep_array(0.0, m);
    array[m] real Ijk = rep_array(0.0, m);

    for (it in 1:max_iter) {
      int idx = 0;
      for (k in 2:p) {
        for (j in 1:(k - 1)) {
          idx += 1;
          int jk = pair_index(j, k, p);
          real pre_j = x[j] - Ij[idx];
          real pre_k = x[k] - Ik[idx];
          real pre_jk = x[jk] - Ijk[idx];
          real post_j = pre_j;
          real post_k = pre_k;
          real post_jk = pre_jk;
          real s = -w[j] * pre_j - w[k] * pre_k + 2 * w[jk] * pre_jk;
          if (s > 0) {
            real t = s / (square(w[j]) + square(w[k]) + 4 * square(w[jk]));
            post_j += t * w[j];
            post_k += t * w[k];
            post_jk -= 2 * t * w[jk];
          }
          Ij[idx] = post_j - pre_j;
          Ik[idx] = post_k - pre_k;
          Ijk[idx] = post_jk - pre_jk;
          x[j] = post_j;
          x[k] = post_k;
          x[jk] = post_jk;
        }
      }
    }
    return x;
  }
}
data {
  int<lower=1> N;
  vector[N] y;
  int<lower=1> p;              // number of main effects
  int<lower=1> K;              // p + p * (p - 1) / 2
  matrix[N, K + 1] X;          // first column is the intercept
  int<lower=1> Ntest;
  vector[Ntest] ytest;
  matrix[Ntest, K + 1] Xtest;
  real<lower=0> R2_mean;
  real<lower=0> R2_prec;
  real<lower=0> proximal_lambda;
  int<lower=1> max_iter;       // Dykstra iterations
  real<lower=0> scale_sigma;
}
transformed data {
  matrix[N, K] Xc;             // centred predictors without the intercept
  vector[K] means_X;
  vector[K] sds_X;             // sd(X_j) for main effects, sqrt(E[(X_j X_k)^2]) for interactions
  vector[K] gw;                // constraint weights 1 / sds_X^2
  real ymean = mean(y);
  vector[N] yc = y - ymean;

  for (i in 1:K) {
    means_X[i] = mean(X[, i + 1]);
    Xc[, i] = X[, i + 1] - means_X[i];
  }
  for (j in 1:p) {
    sds_X[j] = sd(X[, j + 1]);
    gw[j] = 1.0 / square(sds_X[j]);
  }
  for (k in 2:p) {
    for (j in 1:(k - 1)) {
      int jk = pair_index(j, k, p);
      vector[N] x_jk = X[, j + 1] .* X[, k + 1];
      real m2 = variance(x_jk) + square(mean(x_jk));
      sds_X[jk] = sqrt(m2);
      gw[jk] = 1.0 / m2;
    }
  }
}
parameters {
  vector[K] zbeta;
  vector[K] zphi;
  real Intercept;              // intercept of the centred predictors
  real<lower=0, upper=1> R2;
  real<lower=0> sigma;
}
transformed parameters {
  real tau2 = sigma^2 * R2 / (1 - R2);
  vector[K] logit_phi = dykstra(zphi, gw, p, max_iter);
  real dist2 = dot_self(zphi - logit_phi);
  simplex[K] phi = softmax(logit_phi);
  vector[K] beta = zbeta .* sqrt(phi * tau2) ./ sds_X;
}
model {
  target += normal_id_glm_lpdf(yc | Xc, Intercept, beta, sigma);
  target += beta_lpdf(R2 | R2_mean * R2_prec, (1 - R2_mean) * R2_prec);
  target += std_normal_lpdf(zbeta);
  target += std_normal_lpdf(zphi);
  target += -(0.5 / proximal_lambda) * dist2;
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
