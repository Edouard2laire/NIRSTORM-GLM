%% setup path
addpath(genpath('/Users/edelaire1/Documents/software/gramm-master'))
addpath('/Users/edelaire1/Documents/software/brainstorm3'); 
if ~brainstorm('status')
    brainstorm server
end

%% Load data
sFiles = {...
    'sub-01/sub-01_task-tapping_run-01_preproc/results_NIRS_wMEM___smooth=0.6_DWT_j1___2___3___4___5___6___7___8___9__10_____HbO_261024_1309_low.mat', ...
    'sub-01/sub-01_task-tapping_run-01_preproc/results_NIRS_wMEM___smooth=0.6_DWT_j1___2___3___4___5___6___7___8___9__10_____HbR_261024_1309_low.mat'};


sDataHbO = in_bst_results(sFiles{1}, 1);
sDataHbR = in_bst_results(sFiles{2}, 1);

sEvent = in_bst_data(sDataHbO.DataFile, 'Events').Events;
tapping = sEvent(strcmp({sEvent.label}, 'tapping'));

sSubject = bst_get('Subject', 'sub-01');
sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);

atlas_name = 'User scouts'; ROI_name = 'hand';

sAtlas = sCortex.Atlas(strcmp({sCortex.Atlas.Name}, atlas_name));
sRoi = sAtlas.Scouts(strcmp({sAtlas.Scouts.Label}, ROI_name));


time = sDataHbO.Time; 
data = [ mean(sDataHbO.ImageGridAmp(sRoi.Vertices, :), 1); ...
         mean(sDataHbR.ImageGridAmp(sRoi.Vertices, :), 1) ];


%data_downsample = interp1(time, data', time(1:5:end), 'linear')';       % Downsampled Vector
%time = time(1:5:end);
dt = time(2)-time(1);

%% Initialize model

model       = nst_glm_initialize_model(time);

sHRF = readtable('/Users/edelaire1/Documents/software/epilepsy_hrf_lib/HRFbasis_500Hz.csv');

hrf_time = 0:dt:sHRF.time(end); 

hrf_downsample = interp1(sHRF.time, sHRF.pc1, hrf_time, 'linear');       % Downsampled Vector
X_event = nst_make_event_regressors(tapping, hrf_downsample', time);
model.X         = [model.X X_event./max(X_event)];

hrf_downsample = interp1(sHRF.time, sHRF.pc2, hrf_time, 'linear');       % Downsampled Vector
X_event = nst_make_event_regressors(tapping, hrf_downsample', time);
model.X         = [model.X X_event./max(X_event)];

hrf_downsample = interp1(sHRF.time, sHRF.pc3, hrf_time, 'linear');       % Downsampled Vector
X_event = nst_make_event_regressors(tapping, hrf_downsample', time);
model.X         = [model.X X_event./max(X_event)];

tmp = cellfun(@(x) cellfun(@(y) strjoin({x,y}, '-'), {'pc1', 'pc2', 'pc3'} , 'UniformOutput', false),  {tapping.label}, 'UniformOutput', false);
labels = tmp{:};

model.n_roi     = model.n_roi+length(labels);
model.hrf       = hrf_downsample';
model.reg_names = [model.reg_names labels];
% We apply low and high pass filter
model.accept_filter = [model.accept_filter 3*ones(1,length(labels))];

model = nst_glm_add_regressors(model, 'constant');
model = nst_glm_add_regressors(model, "legendre", 4);
model = nst_glm_apply_filter(model, 'FIR_bp', 0.002, 0.1, 0.002);

%% %% Fit model  and plot results ussing nirstormn fit
fitted_model = nst_glm_fit(model, data' , 'OLS_prewhitening');

beta = fitted_model.B;

figure; hold on;
plot(sHRF.time, beta(1, 1) * sHRF.pc1 +  beta(2, 1) * sHRF.pc2  + beta(3, 1) * sHRF.pc3 , 'r');
plot(sHRF.time, beta(1, 2) * sHRF.pc1 +  beta(2, 2) * sHRF.pc2  + beta(3, 2) * sHRF.pc3 , 'b');
yline(0, '--')

%% With stan


pyenv(Version="/Users/edelaire1/miniconda3/envs/cmdstanpy/bin/python",  ExecutionMode="InProcess")

sys = py.importlib.import_module('sys');
logf = py.open('/tmp/stan_live.log', 'w', pyargs('buffering', int32(1)));  % line-buffered
py.setattr(sys, 'stdout', logf);
py.setattr(sys, 'stderr', logf);


cmdstanpy = py.importlib.import_module('cmdstanpy');

stan_model = cmdstanpy.CmdStanModel(pyargs('stan_file', 'glm_model.stan'));
stan_data = py.dict(pyargs( ...
    'n_voxel', int32(size(data,1)), ...
    'n_frame', int32(size(data,2)), ...
    'n_reg',   int32(size(model.X,2)), ...
    'tr', double(dt), ...
    'Y',       py.numpy.array(double(data)), ...
    'X',       py.numpy.array(double(model.X))));

fit = stan_model.sample(pyargs('data', stan_data, 'parallel_chains', int32(4), 'iter_sampling', int32(1000)));

%% Draw results

A = fit.summary();

vals = double(A.to_numpy());
rows = cellfun(@(x) string(x), cell(A.index.tolist()));
cols = cellfun(@(x) string(x), cell(A.columns.tolist()));

T = array2table(vals, ...
    'RowNames', rows, ...
    'VariableNames', matlab.lang.makeValidName(cols));

disp(T)

%% ssimpledraw of the fit

draws = double(fit.stan_variable('beta_orig'));   % numpy array -> MATLAB double

nHrf = 4000;
figure; hold on;
plot(sHRF.time, squeeze(draws(1:nHrf, 1, 1)) * sHRF.pc1' +  squeeze(draws(1:nHrf, 1, 2)) * sHRF.pc2'  + squeeze(draws(1:nHrf, 1, 3)) * sHRF.pc3' , 'r');
plot(sHRF.time, squeeze(draws(1:nHrf, 2, 1)) * sHRF.pc1' +  squeeze(draws(1:nHrf, 2, 2)) * sHRF.pc2'  + squeeze(draws(1:nHrf, 2, 3)) * sHRF.pc3' , 'b');
yline(0, '--')

%%
idx_voxel  = 1;                                  % voxel whose residuals you plotted

beta = double(fit.stan_variable('beta_orig'));
b = squeeze(mean(beta(:, idx_voxel, :), 1));
res = data(1,:)' - model.X * b;

% empirical ACF from the residuals (as before)
[acf, lags] = acf_base(res, 100);


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

%%
idx_voxel = 1;
maxlag    = 100;
Xm        = model.X;                                   % n_frame x n_reg
% dt must be the sampling interval of data_downsample, i.e. the 'tr' passed to Stan

beta   = double(fit.stan_variable('beta_orig'));       % draws x n_voxel x n_reg
y_pred = double(fit.stan_variable('y_pred'));          % draws x n_voxel x n_frame

B  = squeeze(beta(:, idx_voxel, :));                   % draws x n_reg
YP = squeeze(y_pred(:, idx_voxel, :));                 % draws x n_frame

% empirical ACF of the residuals (posterior mean beta)
res = data_downsample(idx_voxel, :)' - Xm * mean(B, 1)';
[acf, lags] = acf_base(res, maxlag);

% simulated noise, per draw: y_pred - X*beta (same draw), then one ACF per draw
noise   = YP - B * Xm';                                % draws x n_frame
noise   = noise(1:10:end, :);                          % thin to save memory
acf_sim = acf_rows(noise, maxlag);                     % draws x (maxlag+1)

lag_s   = lags(:)' * dt;
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




%% Fit GP on the residual

idx_voxel  = 1;                                  % voxel whose residuals you plotted

beta = double(fit.stan_variable('beta_orig'));
b = squeeze(mean(beta(:, idx_voxel, :), 1));
res = data(1,:)' - model.X * b;
res = res(1:400);

stan_model = cmdstanpy.CmdStanModel(pyargs('stan_file', 'stan_gp.stan'));
stan_data = py.dict(pyargs( ...
                            'N', int32(length(res)), ...
                            'dt', (double(0.1)),...
                            'y', py.numpy.array(double(res'))));



fit_gp = stan_model.sample(pyargs('data', stan_data, 'parallel_chains', int32(4), 'iter_sampling', int32(1000)));

%%

v = 1;
Yz  = data(v, :)' / std(data(v, :));
beta = double(fit.stan_variable('beta'));
b = squeeze(mean(beta(:, v, :), 1));
res_new = Yz - model.X * b;                    % X_new includes the Legendre columns

dt = 0.1;
x = res_new - mean(res_new); n = numel(x);
F = fft(x);
k = 1:floor(n/2)-1;
P = abs(F(k+1)).^2 / n;
f = k' / (n * dt);

figure; loglog(f, movmean(P, 5)); hold on
xline(0.1, 'r--'); xlabel('frequency (Hz)'); ylabel('power')


%% Hrfs 

HRFs = squeeze(draws(:, 1, 1)) * sHRF.pc1' +  squeeze(draws(:, 1, 2)) * sHRF.pc2'  + squeeze(draws(:, 1, 3)) * sHRF.pc3';
stats_HBO = zeros(2, size(HRFs,1));

for iHRFs = 1:size(HRFs,1)
   [pks,locs] = findpeaks(HRFs(iHRFs, :), sHRF.time, 'SortStr','descend');
   stats_HBO(:, iHRFs) = [pks(1) ; locs(1)];
end

HRFs = squeeze(draws(:, 2, 1)) * sHRF.pc1' +  squeeze(draws(:, 2, 2)) * sHRF.pc2'  + squeeze(draws(:, 2, 3)) * sHRF.pc3';
stats_HBR = zeros(2, size(HRFs,1));

for iHRFs = 1:size(HRFs,1)
   [pks,locs] = findpeaks(-HRFs(iHRFs, :), sHRF.time, 'SortStr','descend');
   stats_HBR(:, iHRFs) = [-pks(1) ; locs(1)];
end

HRFs = [squeeze(draws(:, 1, 1)) * sHRF.pc1' +  squeeze(draws(:, 1, 2)) * sHRF.pc2'  + squeeze(draws(:, 1, 3)) * sHRF.pc3'; ...
        squeeze(draws(:, 2, 1)) * sHRF.pc1' +  squeeze(draws(:, 2, 2)) * sHRF.pc2'  + squeeze(draws(:, 2, 3)) * sHRF.pc3'];

HRFs_type = [repmat({'HbO'}, 1, 4000), repmat({'HbR'}, 1, 4000)];



%%

figure;
 hold on;
plot(sHRF.time, squeeze(draws(1:10, 1, 1)) * sHRF.pc1' +  squeeze(draws(1:10, 1, 2)) * sHRF.pc2'  + squeeze(draws(1:10, 1, 3)) * sHRF.pc3' , 'r');
plot(sHRF.time, squeeze(draws(1:10, 2, 1)) * sHRF.pc1' +  squeeze(draws(1:10, 2, 2)) * sHRF.pc2'  + squeeze(draws(1:10, 2, 3)) * sHRF.pc3' , 'b');
histogram2(stats_HBO(2, :),stats_HBO(1, :),'DisplayStyle','tile','ShowEmptyBins','on');
histogram2(stats_HBR(2, :),stats_HBR(1, :),'DisplayStyle','tile','ShowEmptyBins','on');
yline(0, '--')


%%

dt_hight = sHRF.time(2) - sHRF.time(1) ;                         % HRF sampling interval (s), must match sHRF's time grid
T  = 35;                          % simulated duration (s)
n_t = round(T / dt_hight);
t   = (0:n_t-1) * dt_hight;

% make the HRF exactly n_t samples (frames past 35 s can't affect the first 35 s)
L = size(HRFs, 2);
if L >= n_t
    hrf = HRFs(:, 1:n_t);
else
    hrf = [HRFs, zeros(size(HRFs,1), n_t - L)];
end

% 10 s boxcar starting at t = 0
stim = double(t < 15);            % 1 x n_t

% causal convolution for all draws at once
S = toeplitz(stim(:), [stim(1), zeros(1, n_t-1)]);   % n_t x n_t
resp = (hrf * S') * dt;           % n_draws x n_t

%% Estimate model prediction

data_pred = zeros(size(draws, 1), 2, size(model.X, 1));
for iS = 1:2
    data_pred(:, iS, :) = squeeze(draws(:, iS, :)) * model.X';
end

y_pred = [squeeze(mean(data_pred(:, 1, :))) , squeeze(mean(data_pred(:, 2, :))), data(1, :)', data(2, :)'];


%%
clear G;


G(1,1) = gramm('x', sHRF.time, 'y', HRFs, 'color', HRFs_type);
G(1,1).stat_summary("geom","area", "type","95percentile");
G(1,1).geom_hline('yintercept',0, 'style', '--');
G(1,1).set_title('HRF time course');

G(1,2) = gramm('x', t, 'y', resp, 'color', HRFs_type);
G(1,2).stat_summary("geom","area", "type","95percentile");
G(1,2).geom_hline('yintercept',0, 'style', '--');
G(1,2).set_title('Response to a 15s tapping block');
G(1,2).axe_property('XLim',[0, 35]);

stats_hrf = [stats_HBO, stats_HBR];


G(2,1) = gramm('x', stats_hrf(2, :)',  'color', HRFs_type);
G(2,1).stat_density();
G(2,1).set_title('Time of peak (s)')


G(2,2) = gramm('x', model.time',  'y', y_pred', 'color', {'HbO - pred'; 'HbR - pred'; 'HbO'; 'HbR'});
G(2,2).geom_line();
G(2,2).geom_hline('yintercept',0, 'style', '--');
G(2,2).axe_property('XLim',[0, model.time(end)]);
G(2,2).set_title('Time course prediction')

figure;
G.draw() 

G(1,1).update('x', stats_hrf(2, :)', 'y', stats_hrf(1, :)', 'color', HRFs_type);
G(1,1).stat_ellipse('type','95percentile','geom','area','patch_opts',{'FaceAlpha',0.1,'LineWidth',2});
G(1,1).axe_property('XLim',[0, 25]);

figure;
G.draw() 
function A = acf_rows(X, maxlag)          % ACF of each row, via FFT
    X = X - mean(X, 2);
    n = size(X, 2);
    nfft = 2^nextpow2(2*n - 1);
    F = fft(X, nfft, 2);
    c = real(ifft(F .* conj(F), nfft, 2));
    A = c(:, 1:maxlag+1) ./ c(:, 1);
end
