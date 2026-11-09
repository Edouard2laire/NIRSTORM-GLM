data {
  int<lower=1> n_voxel;
  int<lower=3> n_frame;
  int<lower=1> n_reg;
  array[n_voxel] vector[n_frame] Y;
  matrix[n_frame, n_reg] X;
}

transformed data {
  vector[n_voxel] y_sd;
  array[n_voxel] vector[n_frame] Y_z;
  for (v in 1:n_voxel) {
    y_sd[v] = sd(Y[v]);
    if (y_sd[v] <= 0) reject("Voxel ", v, " has zero variance");
    Y_z[v] = Y[v] / y_sd[v];
  }
}

parameters {
  array[n_voxel] vector[n_reg] beta;
  array[n_voxel] vector[2] phi;
  vector<lower=0>[n_voxel] sigma_err;
}

model {
  for (v in 1:n_voxel) {
    beta[v] ~ normal(0, 5);
    phi[v] ~ normal(0, 0.5);
  }
  sigma_err ~ normal(0, 10);

  for (v in 1:n_voxel) {
    vector[n_frame] e = Y_z[v] - X * beta[v];
    e[3:n_frame] ~ normal(phi[v][1] * e[2:(n_frame - 1)]
                          + phi[v][2] * e[1:(n_frame - 2)],
                          sigma_err[v]);
  }
}

generated quantities {
  array[n_voxel] vector[n_reg] beta_orig;
  vector[n_voxel] sigma_orig = sigma_err .* y_sd;
  for (v in 1:n_voxel) {
    beta_orig[v] = beta[v] * y_sd[v];
  }
}