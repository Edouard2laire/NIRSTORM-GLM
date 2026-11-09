%% setup path
addpath(genpath('/NAS/home/edelaire/Documents/software/gramm-master'))
addpath('/NAS/home/edelaire/Documents/software/brainstorm3'); 
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
dt = time(2)-time(1);

%% Initialize model
[model, sHRF] = nst_create_model(time, tapping);

%% %% Fit model  and plot results ussing nirstormn fit
fitted_model = nst_glm_fit(model, data' , 'OLS_prewhitening');

beta = fitted_model.B;

figure; hold on;
plot(sHRF.time, beta(1, 1) * sHRF.pc1 +  beta(2, 1) * sHRF.pc2  + beta(3, 1) * sHRF.pc3 , 'r');
plot(sHRF.time, beta(1, 2) * sHRF.pc1 +  beta(2, 2) * sHRF.pc2  + beta(3, 2) * sHRF.pc3 , 'b');
yline(0, '--')

%% With stan

pyenv("Version", '/home/edelaire/anaconda3/envs/cmdstan/bin/python', 'ExecutionMode', 'OutOfProcess')

[stan_model, fit] = fit_model(model, data);

%% Draw results

cellfun(@(x)display_summary(x), fit, UniformOutput=false);

%% Simple draw of the fit

draws_HbO = double(fit{1}.stan_variable('beta_orig'));   % numpy array -> MATLAB double
draws_HbR = double(fit{2}.stan_variable('beta_orig'));   % numpy array -> MATLAB double

nHrf = 30;
figure; hold on;
plot(sHRF.time, squeeze(draws_HbO(1:nHrf, 1)) * sHRF.pc1' +  squeeze(draws_HbO(1:nHrf,  2)) * sHRF.pc2'  + squeeze(draws_HbO(1:nHrf, 3)) * sHRF.pc3' , 'r');
plot(sHRF.time, squeeze(draws_HbR(1:nHrf, 1)) * sHRF.pc1' +  squeeze(draws_HbR(1:nHrf,  2)) * sHRF.pc2'  + squeeze(draws_HbR(1:nHrf, 3)) * sHRF.pc3' , 'b');
yline(0, '--')


%% Hrfs 

HRFs = compute_hrf(sHRF, fit{1});
stats_HBO = zeros(2, size(HRFs,1));

for iHRFs = 1:size(HRFs,1)
   [pks,locs] = max(HRFs(iHRFs, :));
   stats_HBO(:, iHRFs) = [pks(1) ;  sHRF.time(locs(1))];
end

HRFs = compute_hrf(sHRF, fit{2});
stats_HBR = zeros(2, size(HRFs,1));

for iHRFs = 1:size(HRFs,1)
   [pks,locs] = min(HRFs(iHRFs, :));
   stats_HBR(:, iHRFs) = [pks(1) ; sHRF.time(locs(1))];
end

HRFs = [ compute_hrf(sHRF, fit{1}) ;
         compute_hrf(sHRF, fit{2})];
HRFs_type = [repmat({'HbO'}, 1, 4000), repmat({'HbR'}, 1, 4000)];

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

figure;
ax = subplot(2,2,1);
plot_prediction(ax, model, data(1, :), fit{1});
title('Prediction HbO')
xlabel('Time (s)')

ax = subplot(2,2,2);
plot_residual(ax, model, data(1, :), fit{1});
plot_residual(ax, model, data(2, :), fit{2});
title('Residuals'); legend({'HbO', 'CI HbO', 'HbR', 'CI HbR'})
xlabel('Time (s)')

ax = subplot(2,2,3);
plot_prediction(ax, model, data(2, :), fit{2});
title('Prediction HbR'); xlabel('Time (s)')


ax = subplot(2,2,4);
plot_ACF_residual(ax, model, data(1,:), fit{1})
plot_ACF_residual(ax, model, data(2,:), fit{2})
title('Autocorellation function'); legend({'HbO',   'HbR'})
xlabel('Time (s)')

sgtitle('Model evaluation')


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
G(2,1).stat_density("bandwidth", 0.25, 'extra_x',true); 
G(2,1).axe_property('XLim',[7, 9]);
G(2,1).set_title('Time of peak (s)')

G(2,2) = gramm('x', abs(stats_hrf(1, :))',  'color', HRFs_type);
G(2,2).stat_density("bandwidth", 0.25, 'extra_x',true); 
G(2,2).set_title('Value of peak (\muml/l)')
G(2,2).axe_property('XLim',[0, 40]);


figure;
G.draw() 

G(1,1).update('x', stats_hrf(2, :)', 'y', stats_hrf(1, :)', 'color', HRFs_type);
G(1,1).stat_ellipse('type','95percentile','geom','area','patch_opts',{'FaceAlpha',0.1,'LineWidth',2});
G(1,1).axe_property('XLim',[0, 25]);

figure;
G.draw() 
