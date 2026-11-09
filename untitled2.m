pyenv(Version="/Users/edelaire1/miniconda3/envs/cmdstanpy/bin/python",  ExecutionMode="InProcess")


cmdstanpy = py.importlib.import_module('cmdstanpy');
stan_file = fullfile(string(cmdstanpy.cmdstan_path()), 'examples', 'bernoulli', 'bernoulli.stan');

model = cmdstanpy.CmdStanModel(pyargs('stan_file', stan_file));

name = string(model.name)
exe  = string(py.str(model.exe_file))
code = string(model.code())


%% Sample from the model
data = py.dict(pyargs('N', int32(10), 'y', py.list(int32([0 1 0 0 0 0 0 0 0 1]))));

% equivalent Python: model.sample(data=data, chains=4, iter_sampling=2000)
fit = model.sample(pyargs('data', data, 'chains', int32(4), 'iter_sampling', int32(2000)));

draws = double(fit.stan_variable('theta'));   % numpy array -> MATLAB double


figure;
histogram(draws)

%% Display summary 

A = fit.summary();

vals = double(A.to_numpy());
rows = cellfun(@(x) string(x), cell(A.index.tolist()));
cols = cellfun(@(x) string(x), cell(A.columns.tolist()));

T = array2table(vals, ...
    'RowNames', rows, ...
    'VariableNames', matlab.lang.makeValidName(cols));

disp(T)


%% Display diagnoses

disp(string(fit.diagnose()))


%% Save output
 
mkdir('tmp')
fit.save_csvfiles(pyargs('dir', 'tmp'))


%% Draws

draws = double(fit.draws());
names = string(cell(fit.column_names));

idx = find(names == "theta");
theta = draws(:, :, idx);          % 2000 x 4 (draws x chains)

figure;
plot(1:size(theta, 1),  theta)


%% Variational inference

bern_vb = model.variational(pyargs('data', data))

samples = double(bern_vb.variational_sample);


%% 
%%
colnames = string(fit.column_names);
draws = double(fit.draws());
figure;
plot(draws(:, :,  find(colnames == 'beta[1]')))
legend({'Chain #1', 'Chain #2','Chain #3','Chain #4',})

A = draws(:, :,  find(colnames == 'beta[1]'));
B= draws(:, :,  find(colnames == 'beta[2]'));

figure; plot(A(:), B(:), '+')

%%

%%
idx_voxel  = 1;                                  % voxel whose residuals you plotted

beta = double(fit.stan_variable('beta_orig'));
res = data(1,:)' - model.X * mean(beta, 1)';

% empirical ACF from the residuals (as before)
[acf, lags] = acf_base(res, 100);

figure; plot(lags, acf)

lag_s = lags(:)';

y_noise = double(fit.stan_variable('y_noise'));
y_noise = squeeze(y_noise(:, idx_voxel, :));

acf_sim = acf_rows(y_noise, 100);

acf_med = median(acf_sim, 1);
acf_ci  = prctile(acf_sim, [2.5 97.5], 1);

figure; hold on
fill([lag_s fliplr(lag_s)], [acf_ci(1,:) fliplr(acf_ci(2,:))], ...
    [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.6)
plot(lag_s, acf_med, 'r', 'LineWidth', 2)
stem(lag_s, acf, 'filled', 'Color', [0 0.45 0.74], 'MarkerSize', 3)
yline(1.96/sqrt(numel(res)), 'k--'); yline(-1.96/sqrt(numel(res)), 'k--');
xlabel('lag (s)'); ylabel('autocorrelation')
legend('AR(1) simulated, 95% interval', 'AR(1) simulated, median', 'residual ACF', ...
    'Location', 'northeast')

%%

figure; 
subplot(2,1,1)
plot(model.time, y_noise(1:10, :));
title('Simulated AR(1)'); xlabel('Time(s)')
subplot(2,1,2)
plot(model.time, res);
title('Residual'); xlabel('Time(s)')