functions {
  // spectral density of a Matern 5/2 kernel (angular frequency omega)
  vector spd_matern52(vector omega, real alpha, real ell) {
    real a = sqrt(5) / ell;
    vector[rows(omega)] q = square(a) + square(omega);
    return square(alpha) * (16.0 / 3.0) * pow(a, 5) * inv(q .* q .* q);
  }
}

data {
  int<lower=1> n_voxel;
  int<lower=2> n_frame;
  int<lower=1> n_reg;
  int<lower=2> M;                         // number of basis functions, e.g. 120
  real<lower=0> dt;                       // 0.1 s
  real<lower=1> c;                        // boundary factor, e.g. 1.5
  array[n_voxel] vector[n_frame] Y;
  matrix[n_frame, n_reg] X;
}

transformed data {
  real L = c * n_frame * dt / 2;
  vector[n_frame] t = (linspaced_vector(n_frame, 1, n_frame) - (n_frame + 1) / 2.0) * dt;
  matrix[n_frame, M] Phi;
  vector[M] omega;
  for (j in 1:M) {
    omega[j] = pi() * j / (2 * L);
    Phi[:, j] = sqrt(1 / L) * sin(omega[j] * (t + L));
  }

  vector[n_voxel] y_sd;
  matrix[n_frame, n_voxel] Yz;
  for (v in 1:n_voxel) {
    y_sd[v] = sd(Y[v]);
    if (y_sd[v] <= 0) reject("Voxel ", v, " has zero variance");
    Yz[:, v] = Y[v] / y_sd[v];
  }

  // everything that depends on n_frame is computed once here
  matrix[M, M] PtP = crossprod(Phi);
  matrix[M, n_voxel] PtY = Phi' * Yz;
  matrix[M, n_reg] PtX = Phi' * X;
  matrix[n_reg, n_reg] XtX = crossprod(X);
  matrix[n_reg, n_voxel] XtY = X' * Yz;
  vector[n_voxel] YtY;
  for (v in 1:n_voxel) YtY[v] = dot_self(Yz[:, v]);
}

parameters {
  array[n_voxel] vector[n_reg] beta;
  vector<lower=0>[n_voxel] alpha;         // GP marginal sd (z-scale)
  vector<lower=0>[n_voxel] ell;           // GP length scale (s)
  vector<lower=1e-3>[n_voxel] sigma_n;    // white-noise nugget (z-scale)
}

model {
  for (v in 1:n_voxel) beta[v] ~ normal(0, 5);
  alpha ~ normal(0, 2);
  ell ~ lognormal(log(5), 1);
  sigma_n ~ normal(0, 0.2);

  for (v in 1:n_voxel) {
    // residual r = Yz - X*b has covariance sigma^2 I + Phi diag(S) Phi'
    real s2 = square(sigma_n[v]);
    vector[M] sq = sqrt(spd_matern52(omega, alpha[v], ell[v]));
    vector[n_reg] b = beta[v];

    real rr = YtY[v] - 2 * dot_product(b, XtY[:, v]) + quad_form(XtX, b);
    vector[M] u = PtY[:, v] - PtX * b;

    // Woodbury: B = I + D^(1/2) PtP D^(1/2) / s2
    matrix[M, M] LB = cholesky_decompose(add_diag(quad_form_diag(PtP, sq) / s2, 1));
    vector[M] w = mdivide_left_tri_low(LB, sq .* u / sigma_n[v]);

    real quad = (rr - dot_self(w)) / s2;
    real logdet = n_frame * log(s2) + 2 * sum(log(diagonal(LB)));
    target += -0.5 * (quad + logdet);
  }
}

generated quantities {
  array[n_voxel] vector[n_reg] beta_orig;
  vector[n_voxel] alpha_orig = alpha .* y_sd;
  vector[n_voxel] sigma_n_orig = sigma_n .* y_sd;
  for (v in 1:n_voxel) beta_orig[v] = beta[v] * y_sd[v];
}