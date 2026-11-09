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
