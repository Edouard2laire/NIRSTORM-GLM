data {
  int<lower=1> n_voxel;
  int<lower=2> n_frame;
  int<lower=1> n_reg;
  real<lower=0> tr;                       // TR in seconds (0.1)
  array[n_voxel] vector[n_frame] Y;       // raw (demeaned) data
  matrix[n_frame, n_reg] X;
}

transformed data {
  vector[n_voxel] y_sd;
  array[n_voxel] vector[n_frame] Y_z;
  for (v in 1:n_voxel) {
    y_sd[v] = sd(Y[v]);
    if (y_sd[v] <= 0)
      reject("Voxel ", v, " has zero variance");
    Y_z[v] = Y[v] / y_sd[v];
  }
}

parameters {
  array[n_voxel] vector[n_reg] beta;
  vector<lower=0>[n_voxel] sigma_err;
  vector<lower=0>[n_voxel] tau;           // correlation time (s)
}

model {
  for (v in 1:n_voxel) beta[v] ~ normal(0, 5);
  sigma_err ~ normal(0, 10);
  tau ~ lognormal(log(1.5), 0.2);         // ~95% prior mass in 1.0 to 2.2 s

  for (v in 1:n_voxel) {
    real rho = exp(-tr / tau[v]);
    vector[n_frame] e = Y_z[v] - X * beta[v];

    e[1] ~ normal(0, sigma_err[v] / sqrt(1 - square(rho)));
    e[2:n_frame] ~ normal(rho * e[1:(n_frame - 1)], sigma_err[v]);
  }
}

generated quantities {
  array[n_voxel] vector[n_reg] beta_orig;
  vector[n_voxel] sigma_orig = sigma_err .* y_sd;
  vector[n_voxel] lambda = inv(tau);      // decay rate (1/s)
  array[n_voxel] vector[n_frame] y_noise;  // simulated time course, original data scale

  for (v in 1:n_voxel) {
    real rho = exp(-tr / tau[v]);
    beta_orig[v] = beta[v] * y_sd[v];

    real e_prev = normal_rng(0, sigma_err[v] / sqrt(1 - square(rho)));
    y_noise[v][1] = e_prev;

    for (t in 2:n_frame) {
      real e_t = normal_rng(rho * e_prev, sigma_err[v]);   // AR(1) innovation
      y_noise[v][t] = e_t;
      e_prev = e_t;
    }

    y_noise[v] = y_noise[v] * y_sd[v];      // back to the original scale
  }
}