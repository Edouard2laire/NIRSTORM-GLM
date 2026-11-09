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
  array[n_voxel] vector[n_reg] beta;      // betas on the z-scored scale
  vector<lower=0>[n_voxel] sigma_err;     // innovation sd on the z-scored scale
  vector<lower=0>[n_voxel] lambda;        // decay rate (1/s)
}

model {
  for (v in 1:n_voxel) {
    beta[v] ~ normal(0, 5);
  }
  
  sigma_err ~ normal(0, 10);
  lambda ~ lognormal(0, 1);

  for (v in 1:n_voxel) {
    real rho = exp(-lambda[v] * tr);
    vector[n_frame] e = Y_z[v] - X * beta[v];

    e[1] ~ normal(0, sigma_err[v] / sqrt(1 - square(rho)));
    e[2:n_frame] ~ normal(rho * e[1:(n_frame - 1)], sigma_err[v]);
  }
}

generated quantities {
  // quantities back on the original data scale
  array[n_voxel] vector[n_reg] beta_orig;
  vector[n_voxel] sigma_orig;
  vector[n_voxel] tau = inv(lambda);      // correlation time (s), scale-free

  for (v in 1:n_voxel) {
    beta_orig[v] = beta[v] * y_sd[v];
  }
  sigma_orig = sigma_err .* y_sd;
}