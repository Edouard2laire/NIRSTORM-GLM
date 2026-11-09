functions {
  // Durbin-Levinson: partial autocorrelations in (-1,1) -> stationary AR coefficients
  vector pacf_to_ar(vector r) {
    int p = rows(r);
    vector[p] phi = rep_vector(0, p);
    phi[1] = r[1];
    for (k in 2:p) {
      vector[k - 1] prev = phi[1:(k - 1)];
      phi[1:(k - 1)] = prev - r[k] * reverse(prev);
      phi[k] = r[k];
    }
    return phi;
  }
}

data {
  int<lower=1> n_voxel;
  int<lower=3> n_frame;
  int<lower=1> n_reg;
  int<lower=1, upper=n_frame - 1> p;            // AR order, a data variable
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
  array[n_voxel] vector<lower=-1, upper=1>[p] pacf;
  vector<lower=0>[n_voxel] sigma_err;
}

transformed parameters {
  array[n_voxel] vector[p] phi;
  for (v in 1:n_voxel) phi[v] = pacf_to_ar(pacf[v]);
}

model {
  for (v in 1:n_voxel) {
    beta[v] ~ normal(0, 5);
    pacf[v] ~ normal(0, 0.5);             // weak; the data dominate (lag 1 will be near 1)
  }
  sigma_err ~ normal(0, 1);

  for (v in 1:n_voxel) {
    vector[n_frame] e = Y_z[v] - X * beta[v];
    vector[n_frame - p] mu = rep_vector(0, n_frame - p);
    for (k in 1:p) mu += phi[v][k] * e[(p + 1 - k):(n_frame - k)];
    e[(p + 1):n_frame] ~ normal(mu, sigma_err[v]);   // conditional on the first p frames
  }
}

generated quantities {
  array[n_voxel] vector[n_reg] beta_orig;
  vector[n_voxel] sigma_orig = sigma_err .* y_sd;
  for (v in 1:n_voxel) beta_orig[v] = beta[v] * y_sd[v];
}