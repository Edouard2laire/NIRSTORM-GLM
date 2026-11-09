function plot_ACF_residual(ax, model, data, fit)
    
    maxlag = 100;
    betas = double(fit.stan_variable('beta_orig'));
    residual = data' - model.X * betas';
    
    dt = model.time(2) - model.time(1);

    lags = dt* (0:maxlag);
    A = acf_rows(residual', maxlag);

    hold(ax, 'on')
    plot(ax, lags, mean(A));
    hold(ax, 'off')
end