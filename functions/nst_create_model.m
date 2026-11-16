function [model, sHRF] = nst_create_model(time, sEvents, sHRF)
%NST_CREATE_MODEL initialize model based on sEvents

    model       = nst_glm_initialize_model(time);
    
    if nargin < 3
        sHRF = nst_load_hrf_lib('HRFbasis_500Hz', time);
    end

    for iEvent = 1:length(sEvents)

        for iComponnent = 1:size(sHRF.basis_functions, 2)
            X_event = nst_make_event_regressors(sEvents(iEvent), sHRF.basis_functions(:, iComponnent), time);
            model.X         = [model.X X_event./max(X_event)];


            model.reg_names{end+1} = sprintf('%s - pc%d', sEvents(iEvent).label, iComponnent); 
            model.accept_filter(end+1) = 3;
            model.n_roi     = model.n_roi + 1;
        end
        model.hrf       = sHRF.basis_functions;
    end

    model = nst_glm_add_regressors(model, 'constant');
    model = nst_glm_add_regressors(model, "legendre", 4);

    %model = nst_glm_add_regressors(model, "DCT", [0.2, 0.6], {'respiration'});
    %model = nst_glm_add_regressors(model, "DCT", [0.8, 1.5], {'cardiac'});

end
