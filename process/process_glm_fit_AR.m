function varargout = process_glm_fit_AR( varargin )
% process_glm_fit_AR: Fit GLM model to obtain the HRF of epileptic activity
% during different sleep stage. The HRF is modeled as a linear combination
% of 3 PCA components. Model is fit using Huppert toolbox (assuming AR(p)
% noise). Results for the subject are obtained using a fixed effect model
% 
% 
%
% Authors: Edouard Delaire, 2026

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() 
    % Description the process
    sProcess.Comment     = 'GLM - cpmpute GLM';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'EPI', 'GLM'};
    sProcess.Index       = 2300;
    sProcess.isSeparator = 0;
    sProcess.Description = '';

    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'results'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.nOutputs    = 1;
    
    sProcess.options.label3.Comment = '<U><B>Design Matrix</B></U>:';
    sProcess.options.label3.Type    = 'label';
    
    sProcess.options.stim_events.Comment = 'Stimulation events: ';
    sProcess.options.stim_events.Type    = 'text';
    sProcess.options.stim_events.Value   = '';

    sProcess.options.out_name.Comment = 'Output Folder name';
    sProcess.options.out_name.Type='text';
    sProcess.options.out_name.Value = '';

end



%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment; 
end

function OutputFiles = Run(sProcess, sInputs) 

    OutputFiles = {};

    %% Select events
    if isempty(sProcess.options.stim_events.Value)
         bst_error('No event selected');
         return;
    end
    selected_event_names = cellfun(@strtrim, strsplit(sProcess.options.stim_events.Value, ','),...
                                   'UniformOutput', 0);

    % Load cortex
    subjectName = sInputs(1).SubjectName;
    sSubject = bst_get('Subject', sInputs(1).SubjectName);
    sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);

    atlas_name = 'User scouts'; ROI_name = 'ReconFocal';

    sAtlas = sCortex.Atlas(strcmp({sCortex.Atlas.Name}, atlas_name));
    sRoi = sAtlas.Scouts(strcmp({sAtlas.Scouts.Label}, ROI_name));

    sFiles = {sInputs.FileName};
    data_types = {'HbO','HbR'};

    iStudy = db_add_condition(subjectName, sProcess.options.out_name.Value);
    sStudy = bst_get('Study', iStudy);

    for iChromophore = 1:length(data_types)

        data_type = data_types{iChromophore};
        dataFiles = sFiles(contains(lower(sFiles), lower(data_type)));

        bst_progress('text', sprintf('Fitting the model for %s', data_type)); 
        [sHRF, models, sDataAll] = initialize_models(dataFiles, selected_event_names);


        sDataOut = repmat(db_template('resultsmat'), 1, length(selected_event_names));
        for iEvent = 1:length(selected_event_names)
            sDataOut(iEvent).ImageGridAmp = zeros(size(sCortex.Vertices, 1), length(sHRF.time));
            sDataOut(iEvent).Std = zeros(size(sCortex.Vertices, 1), length(sHRF.time));
            sDataOut(iEvent).Time = sHRF.time;
            sDataOut(iEvent).Comment = sprintf('HRF %s - %s', selected_event_names{iEvent}, data_type);
            sDataOut(iEvent).DisplayUnits = 'mumol.l-1';
            sDataOut(iEvent).SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
        end


        sDataOut_beta = db_template('resultsmat');
        sDataOut_beta.ImageGridAmp = zeros(size(sCortex.Vertices, 1), length(selected_event_names)*sHRF.nComponnent);
        sDataOut_beta.Std = zeros(size(sCortex.Vertices, 1), length(selected_event_names)*sHRF.nComponnent, length(selected_event_names)*sHRF.nComponnent);
        sDataOut_beta.Time = 1:(length(selected_event_names)*sHRF.nComponnent);
        sDataOut_beta.Comment = sprintf('beta %s - %s', selected_event_names{iEvent}, data_type);
        sDataOut_beta.stats = [];
        sDataOut_beta.DisplayUnits = 'mumol.l-1';
        sDataOut_beta.SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;


        stats = [];
        for iData = 1:length(sDataAll)
            data = getDataRoi(sDataAll(iData), sRoi);
            stats = [stats,nirs.math.ar_irls(data', models(iData).X, 10, 4.685, 0, 0, 1 )];
        end

        stats_fixed = nst_glm_fixed_effect(stats, models);
        stats_fixed.ROI = sRoi;
        HRFs = compute_HRF(stats_fixed, sHRF, selected_event_names);

        for iEvent = 1:length(selected_event_names)
            sDataOut(iEvent).ImageGridAmp(sRoi.Vertices, :) = repmat(mean(HRFs(iEvent).draw), length(sRoi.Vertices), 1);
            sDataOut(iEvent).Std(sRoi.Vertices, :) = repmat(std(HRFs(iEvent).draw), length(sRoi.Vertices), 1);
            sDataOut(iEvent).Leff  =  HRFs(iEvent).nEff;
            sDataOut(iEvent).nAvg   =  HRFs(iEvent).nEff;
        end

        % Generate a new file name in the same folder
        for iEvent = 1:length(selected_event_names)
            OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), sprintf('results_hrf_%s_%s', nst_protect_fn_str(selected_event_names{iEvent}), data_type));
            bst_save(OutputFile, sDataOut(iEvent), 'v7');
            % Register in database
            db_add_data(iStudy, OutputFile, sDataOut(iEvent));
            OutputFiles{end+1} = OutputFile;
        end

        sDataOut_beta.ImageGridAmp(sRoi.Vertices, :) = repmat(stats_fixed.beta', length(sRoi.Vertices), 1);
        sDataOut_beta.Std(sRoi.Vertices, :, :) = repmat(reshape(stats_fixed.covb, [1, size(stats_fixed.covb)]), length(sRoi.Vertices), 1, 1);
        sDataOut_beta.stats = [sDataOut_beta.stats, stats_fixed];
        sDataOut_beta.Leff  =  stats_fixed.nEff;
        sDataOut_beta.nAvg   =  stats_fixed.nEff;
        
        OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), sprintf('results_stats_%s_%s', nst_protect_fn_str(selected_event_names{iEvent}), data_type));
        bst_save(OutputFile, sDataOut_beta, 'v7');
        % Register in database
        db_add_data(iStudy, OutputFile, sDataOut_beta);
        OutputFiles{end+1} = OutputFile;

    end

end

function [sHRF, models, sDataAll] = initialize_models(sFiles, eventNames)
    
    sHRF = nst_load_hrf_lib('HRFbasis_500Hz', 0.1);
    models = []; sDataAll = [];

    for iData = 1:length(sFiles)
        sData = in_bst_results(sFiles{iData}, 0, 'Time', 'DataFile', 'ImageGridAmp');    
        sEvents = in_bst_data(sData.DataFile, 'Events').Events;

        [events, missing_events] = GetEvents(sEvents, eventNames);
        [model, sHRF] = nst_create_model(sData.Time, events, sHRF);
        model = nst_glm_apply_filter(model, 'FIR_bp', 0.002, 0.1, 0.002);

        % save splitted data
        models = [models, model];
        sDataAll = [sDataAll, sData];
    end
end        

function data = getDataRoi(sData, sRoi)

    F = sData.ImageGridAmp;

    if ~iscell(F)
        F = {F};
    end
    
    vertex = sRoi.Vertices;

    % create spatial mapping and reduce the spatial dimension 
    mapping = sparse(1:length(vertex), vertex, 1, length(vertex), size(F{1}, 1));
    F{1}  = mapping * F{1};


    data = bst_multiply_cellmat(F);
    if length(vertex) > 1
        data = mean(data);
    end

end

function HRFs = compute_HRF(stats, sHRF, eventNames)

    draws = mvnrnd(stats.beta, stats.covb, 2000); 
    basis_function = sHRF(1).basis_functions';

    HRFs = repmat(struct( 'eventName', '', 'time', [], 'draw', [], 'nEff', 0), 1, length(eventNames));

    for iEvent = 1:length(eventNames)
        iReg = contains(stats.reg_names, eventNames(iEvent));
    
        if ~any(iReg)
            continue;
        end

        HRFs(iEvent).time = sHRF.time; 
        HRFs(iEvent).draw = draws(:, iReg)*basis_function; 
        HRFs(iEvent).eventName = eventNames{iEvent}; 
        HRFs(iEvent).nEff = unique(stats.nEff(iReg)); 
    end

end