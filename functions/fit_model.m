function [stan_model, fit] = fit_model(model, data)
%FIT_MODEL Summary of this function goes here
%   Detailed explanation goes here
    
    cmdstanpy = py.importlib.import_module('cmdstanpy');
    stan_model = cmdstanpy.CmdStanModel(pyargs('stan_file', fullfile('models','glm_model.stan')));
    
    fit = cell(1, size(data, 1));
    for iData = 1:size(data, 1)
        stan_data = py.dict(pyargs( ...
            'n_frame', int32(size(data,2)), ...
            'n_reg',   int32(size(model.X,2)), ...
            'Y',       py.list(double(data(iData, :))), ...
            'X',       py.numpy.array(double(model.X))));
        
        fit{iData} = stan_model.sample(pyargs('data', stan_data, 'show_console','true',  'parallel_chains', int32(4), 'threads_per_chain', int32(5), 'max_treedepth', int32(20), 'iter_sampling', int32(1000)));
    end
    
end