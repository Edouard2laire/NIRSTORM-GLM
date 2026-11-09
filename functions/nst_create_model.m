function [model, sHRF] = nst_create_model(time, sEvents)
%NST_CREATE_MODEL initialize model based on sEvents

    model       = nst_glm_initialize_model(time);
    
    hrf_file = fullfile(fileparts(which("nst_create_model.m")), '..', 'lib', 'HRFbasis_500Hz.csv');
    if ~file_exist(hrf_file)
        error('Unable to find the HrF library')
    end

    sHRF = readtable(hrf_file);

    dt= time(2) - time(1);
    hrf_time = 0:dt:sHRF.time(end); 

    % Downsample hrf
    pc1 = interp1(sHRF.time, sHRF.pc1, hrf_time, 'linear')';       % Downsampled Vector
    pc2 = interp1(sHRF.time, sHRF.pc2, hrf_time, 'linear')';       % Downsampled Vector
    pc3 = interp1(sHRF.time, sHRF.pc3, hrf_time, 'linear')';       % Downsampled Vector

    sHRF = table(hrf_time', pc1, pc2, pc3, 'VariableNames',{'time', 'pc1', 'pc2', 'pc3'});

    for iEvent = 1:length(sEvents)
        X_event = nst_make_event_regressors(sEvents, sHRF.pc1, time);
        model.X         = [model.X X_event./max(X_event)];
        
        X_event = nst_make_event_regressors(sEvents, sHRF.pc2, time);
        model.X         = [model.X X_event./max(X_event)];
        
        X_event = nst_make_event_regressors(sEvents,  sHRF.pc3, time);
        model.X         = [model.X X_event./max(X_event)];
        
        tmp = cellfun(@(x) cellfun(@(y) strjoin({x,y}, '-'), {'pc1', 'pc2', 'pc3'} , 'UniformOutput', false),  {sEvents.label}, 'UniformOutput', false);
        labels = tmp{:};
        
        model.n_roi     = model.n_roi+length(labels);
        model.hrf       = pc1;
        model.reg_names = [model.reg_names labels];
        % We apply low and high pass filter
        model.accept_filter = [model.accept_filter 3*ones(1,length(labels))];
    end

    model = nst_glm_add_regressors(model, 'constant');
    model = nst_glm_add_regressors(model, "legendre", 4);
    model = nst_glm_apply_filter(model, 'FIR_bp', 0.002, 0.1, 0.002);
end
