%% setup path
addpath(genpath('/NAS/home/edelaire/Documents/software/gramm-master'))
addpath('/NAS/home/edelaire/Documents/software/brainstorm3'); 
if ~brainstorm('status')
    brainstorm server
end

% Load cortex
sSubject = bst_get('Subject', 'PA07');
sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);

atlas_name = 'User scouts'; ROI_name = 'ReconFocal';

sAtlas = sCortex.Atlas(strcmp({sCortex.Atlas.Name}, atlas_name));
sRoi = sAtlas.Scouts(strcmp({sAtlas.Scouts.Label}, ROI_name));


%% Specify data
sFiles{1} = {...
    'PA07/sub-PA07_ses-EEGfNIRS03_run-02_task-sleep_mod-nirs_preproc/results_NIRS_wMEM___smooth=0.6_DWT_j1___2___3___4___5___6___7___8___9__10__11_____HbO___nbo_=_6_selected_261103_1751.mat', ...
    'PA07/sub-PA07_ses-EEGfNIRS03_run-02_task-sleep_mod-nirs_preproc/results_NIRS_wMEM___smooth=0.6_DWT_j1___2___3___4___5___6___7___8___9__10__11_____HbR___nbo_=_6_selected_261103_1751.mat'};
sFiles{2} = {...
    'PA07/sub-PA07_ses-EEGfNIRS03_run-02_task-sleep_mod-nirs_preproc/results_NIRS_wMEM___smooth=0.6_DWT_j1___2___3___4___5___6___7___8___9__10__11_____HbO___nbo_=_6_selected_261116_1606.mat', ...
    'PA07/sub-PA07_ses-EEGfNIRS03_run-02_task-sleep_mod-nirs_preproc/results_NIRS_wMEM___smooth=0.6_DWT_j1___2___3___4___5___6___7___8___9__10__11_____HbR___nbo_=_6_selected_261116_1606.mat'};
sFiles{3} = {...
    'PA07/sub-PA07_ses-EEGfNIRS03_run-02_task-sleep_mod-nirs_preproc/results_NIRS_wMEM___smooth=0.6_DWT_j1___2___3___4___5___6___7___8___9__10__11__12_____HbO_selected_261116_1607.mat', ...
    'PA07/sub-PA07_ses-EEGfNIRS03_run-02_task-sleep_mod-nirs_preproc/results_NIRS_wMEM___smooth=0.6_DWT_j1___2___3___4___5___6___7___8___9__10__11__12_____HbR_selected_261116_1607.mat'};
sFiles{4} = {...
    'PA07/sub-PA07_ses-EEGfNIRS03_run-02_task-sleep_mod-nirs_preproc/results_NIRS_wMEM___smooth=0.6_DWT_j1___2___3___4___5___6___7___8___9__10__11__12_____HbO_selected_261116_1608.mat', ...
    'PA07/sub-PA07_ses-EEGfNIRS03_run-02_task-sleep_mod-nirs_preproc/results_NIRS_wMEM___smooth=0.6_DWT_j1___2___3___4___5___6___7___8___9__10__11__12_____HbR_selected_261116_1608.mat'};

eventNames = {'spikes^N2', 'spikes^N3',  'bursts^N2_F', 'bursts^N3_F'};

%% Initialize model

sHRF = nst_load_hrf_lib('HRFbasis_500Hz', 0.1);
modelAll = []; models = [];
dataAll = []; datas = {};
for iData = 1:length(sFiles)
    
    sFile = sFiles{iData};

    sDataHbO = in_bst_results(sFile{1}, 1);
    sDataHbR = in_bst_results(sFile{2}, 1);

    time = sDataHbO.Time; 
    data = [    mean(sDataHbO.ImageGridAmp(sRoi.Vertices, :), 1); ...
                mean(sDataHbR.ImageGridAmp(sRoi.Vertices, :), 1) ];
    dt = time(2)-time(1);

    sEvents = in_bst_data(sDataHbO.DataFile, 'Events').Events;
    [events, missing_events] = GetEvents(sEvents, eventNames);

    [model, sHRF] = nst_create_model(time, events, sHRF);
    model = nst_glm_apply_filter(model, 'FIR_bp', 0.002, 0.1, 0.002);

    % save concatanated data
    modelAll = nst_glm_concatenate(modelAll, model);
    dataAll = [ data, dataAll];

    % save splitted data
    models = [models, model];
    datas{iData} =  data;

end

%% %% Fit model  and plot results ussing nirstormn fit
fitted_model = nst_glm_fit(model, data(1,:)' , 'OLS_prewhitening');
beta = fitted_model.B;

for iEvent = 1:length(eventNames)
    iReg = contains(model.reg_names, eventNames(iEvent));
    
    if ~any(iReg)
        continue;
    end

    figure;
    hold on;
    plot(sHRF.time, model.hrf * beta(iReg, 1) , 'r');
    %plot(sHRF.time, model.hrf * beta(iReg, 2) , 'b');
    title(eventNames(iEvent))
    yline(0, '--')
end



%% Fit using Huppert toolbox

addpath(genpath('/NAS/home/edelaire/Documents/software/nirs-toolbox'));

stats(1) = nirs.math.ar_irls(dataAll(1,:)', modelAll.X, 10, 4.685, 0, 0, 1 );
stats(2) = nirs.math.ar_irls(dataAll(2,:)', modelAll.X, 10, 4.685, 0, 0, 1 );


draws_HbO = mvnrnd(stats(1).beta, stats(1).covb, 2000); 
draws_HbR = mvnrnd(stats(2).beta, stats(2).covb, 2000); 


for iEvent = 1:length(eventNames)
    iReg = contains(modelAll.reg_names, eventNames(iEvent));

    if ~any(iReg)
        continue;
    end

    y_HbO = draws_HbO(:, iReg)*(sHRF.basis_functions)'; 
    y_HbR = draws_HbR(:, iReg)*(sHRF.basis_functions)'; 

    figure;
    hold on;
    shadedErrorBar(sHRF.time, mean(y_HbO,1), std(y_HbO), 'lineProps', 'r-');
    shadedErrorBar(sHRF.time, mean(y_HbR,1), std(y_HbR), 'lineProps', 'b-');
    title(strrep(eventNames(iEvent), '^', ' '))
    yline(0, '--')
end

%% Fit one model per data

for iData = 1:length(datas)

    data = datas{iData};

    stats_HbO(iData) = nirs.math.ar_irls(data(1,:)', models(iData).X, 10, 4.685, 0, 0, 1 );
    stats_HbR(iData) = nirs.math.ar_irls(data(2,:)', models(iData).X, 10, 4.685, 0, 0, 1 );

end

stats_HbO_all = nst_glm_fixed_effect(stats_HbO, models);
stats_HbR_all = nst_glm_fixed_effect(stats_HbR, models);


draws_HbO = mvnrnd(stats_HbO_all.beta, stats_HbO_all.covb, 2000); 
draws_HbR = mvnrnd(stats_HbR_all.beta, stats_HbR_all.covb, 2000); 


for iEvent = 1:length(eventNames)
    iReg = contains(stats_HbO_all.reg_names, eventNames(iEvent));

    if ~any(iReg)
        continue;
    end

    y_HbO = draws_HbO(:, iReg)*(sHRF.basis_functions)'; 
    y_HbR = draws_HbR(:, iReg)*(sHRF.basis_functions)'; 

    figure;
    hold on;
    shadedErrorBar(sHRF.time, mean(y_HbO,1), std(y_HbO), 'lineProps', 'r-');
    shadedErrorBar(sHRF.time, mean(y_HbR,1), std(y_HbR), 'lineProps', 'b-');
    title(eventNames(iEvent))
    yline(0, '--')
end