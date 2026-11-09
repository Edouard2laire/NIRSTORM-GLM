function HRFs = compute_hrf(sHRF, fit)

    draws_HbO = double(fit.stan_variable('beta_orig'));   % numpy array -> MATLAB double
    HRFs = squeeze(draws_HbO(:, 1)) * sHRF.pc1' +  squeeze(draws_HbO(:, 2)) * sHRF.pc2'  + squeeze(draws_HbO(:,3)) * sHRF.pc3';

end