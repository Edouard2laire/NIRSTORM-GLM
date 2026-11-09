functions {
  // log-likelihood of y = f + noise, f ~ GP(0, Matern 3/2), regular grid spacing dt
  real gp_matern32_ss_lpdf(vector y, real dt, real alpha, real rho, real sigma) {
    int N = num_elements(y);
    real lam = sqrt(3) / rho;
    real e = exp(-lam * dt);
    matrix[2, 2] A = e * [[1 + lam * dt, dt],
                          [-square(lam) * dt, 1 - lam * dt]];
    matrix[2, 2] Pinf = diag_matrix([square(alpha), square(lam * alpha)]');
    matrix[2, 2] Q = Pinf - A * Pinf * A';
    real sq_sigma = square(sigma);

    vector[2] m = rep_vector(0, 2);
    matrix[2, 2] P = Pinf;          // stationary prior at t_1
    real lp = 0;

    for (n in 1:N) {
      // update with observation y[n]
      real S = P[1, 1] + sq_sigma;
      real v = y[n] - m[1];
      vector[2] K = P[:, 1] / S;
      lp += -0.5 * (log(2 * pi() * S) + square(v) / S);
      m = m + K * v;
      P = P - K * P[1, :];
      // predict to next step
      m = A * m;
      P = A * P * A' + Q;
    }
    return lp;
  }
}
data {
  int<lower=1> N;
  real<lower=0> dt;      // 0.1
  vector[N] y;           // standardize before passing in
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}
model {
  rho   ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  sigma ~ std_normal();
  y ~ gp_matern32_ss(dt, alpha, rho, sigma);
}