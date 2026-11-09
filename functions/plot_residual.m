function plot_residual(ax, model, data, fit)

hold(ax, 'on')
betas = mean(double(fit.stan_variable('beta_orig')));
prediction = data' - model.X * betas';

err = mean(double(fit.stan_variable('sigma_orig')));
p = plot(ax, model.time, prediction);

X_fill = [model.time, fliplr(model.time)];
Y_fill = [prediction' - err, fliplr(prediction' + err)];

fill(X_fill, Y_fill, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', p.Color);
xlim(ax, [model.time(1), model.time(end)]);
hold(ax, 'off')

end