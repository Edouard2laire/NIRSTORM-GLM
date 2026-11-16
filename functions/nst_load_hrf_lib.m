function sHRF = nst_load_hrf_lib(name, dt)

    hrf_file = fullfile(fileparts(which("nst_create_model.m")), '..', 'lib', [name, '.csv']);
    if ~file_exist(hrf_file)
        error('Unable to find the HrF library')
    end
    
    % If user provided time vector, compute dt
    if ~isscalar(dt)
        dt = dt(2) - dt(1);
    end

    sHRF = readtable(hrf_file);
    hrf_time = 0:dt:sHRF.time(end); 

    colnames = sHRF.Properties.VariableNames(2:end);
    nComponnent = length(colnames);
    basis_functions = zeros(length(hrf_time), nComponnent);


    % Downsample hrf
    for iComponent = 1:length(colnames)
        basis_functions(:, iComponent) = interp1(sHRF.time, sHRF.(colnames{iComponent}), hrf_time, 'linear')';
    end
    

    sHRF = struct('time', hrf_time, 'basis_functions', basis_functions, 'nComponnent', nComponnent);

end