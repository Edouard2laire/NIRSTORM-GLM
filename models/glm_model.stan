data {
  int<lower=2> n_frame;
  int<lower=1> n_reg;
  vector[n_frame] Y;       // raw (demeaned) data
  matrix[n_frame, n_reg] X;
}

transformed data {

  real y_sd = sd(Y);
  vector[n_frame] Y_z = Y / y_sd;
}

parameters {
  vector[n_reg] beta;
  real<lower=0> sigma_err;
}

model {

  beta ~ normal(0, 5);
  sigma_err ~ normal(0, 10);
  Y_z ~ normal(X * beta, sigma_err);

}

generated quantities {
  vector[n_reg] beta_orig =  beta * y_sd;
  real sigma_orig = sigma_err * y_sd;
}