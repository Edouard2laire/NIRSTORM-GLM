function plot_prediction(ax, model, data, fit)

hold(ax, 'on')
plot(ax, model.time, data);

betas = mean(double(fit.stan_variable('beta_orig')));
prediction = model.X * betas';

err = mean(double(fit.stan_variable('sigma_orig')));


plot(ax, model.time, prediction);

X_fill = [model.time, fliplr(model.time)];
Y_fill = [prediction' - err, fliplr(prediction' + err)];

fill(X_fill, Y_fill, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

xlim(ax, [model.time(1), model.time(end)]);
legend({'Data', 'Prediction'})
hold(ax, 'off')

end