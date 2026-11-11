%% setup path
addpath(genpath('/NAS/home/edelaire/Documents/software/gramm-master'))
addpath('/NAS/home/edelaire/Documents/software/brainstorm3'); 
if ~brainstorm('status')
    brainstorm server
end

%% Load data

sFiles = {...
    'sub-01/sub-01_task-tapping_run-01_preproc/results_NIRS_wMEM___smooth=0.6_DWT_j1___2___3___4___5___6___7___8___9__10_____HbO_261024_1309.mat', ...
    'sub-01/sub-01_task-tapping_run-01_preproc/results_NIRS_wMEM___smooth=0.6_DWT_j1___2___3___4___5___6___7___8___9__10_____HbR_261024_1309.mat'};


sDataHbO = in_bst_results(sFiles{1}, 1);
sDataHbR = in_bst_results(sFiles{2}, 1);

sEvent = in_bst_data(sDataHbO.DataFile, 'Events').Events;
tapping = sEvent(strcmp({sEvent.label}, 'tapping'));

sSubject = bst_get('Subject', 'sub-01');
sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);

atlas_name = 'Atlas #9';

time = sDataHbO.Time; 
dt = time(2)-time(1);

%% Initialize model
[model, sHRF] = nst_create_model(time, tapping);


%%

sAtlas = sCortex.Atlas(strcmp({sCortex.Atlas.Name}, atlas_name));
sRoi = sAtlas.Scouts;

sOutMu = zeros( size(sCortex.Vertices, 1), length(sHRF.time));
sOutSTD = zeros( size(sCortex.Vertices, 1), length(sHRF.time));

for iROI = 1:length(sRoi)
    data = mean(sDataHbR.ImageGridAmp(sRoi(iROI).Vertices, :), 1);

    stats = nirs.math.ar_irls(data', model.X, 10, 4.685, 0, 0, 1 );
    try
        draws = mvnrnd(stats.beta, stats.covb, 2000); 
        y = draws(:,1:3)*([sHRF.pc1 , sHRF.pc2 ,  sHRF.pc3])'; [ValMaxHbO, ImaxHbO] = max(y, [], 2);
    
        sOutMu(sRoi(iROI).Vertices, :) = repmat(mean(y), length(sRoi(iROI).Vertices), 1);
        sOutSTD(sRoi(iROI).Vertices, :) = repmat(std(y), length(sRoi(iROI).Vertices), 1);
    catch
    end
end

sOutBst = db_template('resultsmat');
sOutBst.ImageGridAmp = sOutMu;
sOutBst.Std = sOutSTD/2;
sOutBst.Time = sHRF.time;
sOutBst.Comment = 'GLM result | HbR';
sOutBst.SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;


