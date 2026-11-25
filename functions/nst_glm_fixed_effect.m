function stats_fe = nst_glm_fixed_effect(stats, models, nuisance_patterns)
% NST_GLM_FIXED_EFFECT  Fixed-effect combination of per-run GLM stats.
%
%   stats_fe = nst_glm_fixed_effect(stats, models)
%   stats_fe = nst_glm_fixed_effect(stats, models, nuisance_patterns)
%
%   stats   : struct array (1xN) or cell array (1xN) of outputs of
%             nirs.math.ar_irls, with fields
%               beta (nReg x nChan), covb (nReg x nReg x nChan), dfe
%   models  : struct array / cell array of models (same length as stats),
%             each with a .reg_names field (one name per column of X).
%             Alternatively, a cell array of cell arrays of names.
%   nuisance_patterns : cellstr of substrings identifying nuisance
%             regressors (default: {'Constant', 'Poly'})
%
%   Information regressors are matched by name across runs (zero
%   precision where a run lacks a regressor). Nuisance regressors are
%   treated as per-run and are marginalised out (dropped).
%
%   stats_fe fields: beta, covb, dfe, tstat, pval, reg_names

    if nargin < 3 || isempty(nuisance_patterns)
        nuisance_patterns = {'Constant', 'Poly'};
    end

    % --- Normalise inputs to cell arrays ---
    if isstruct(stats),  stats  = num2cell(stats);  end
    nRun = numel(stats);

    allNames = cell(1, nRun);
    for iRun = 1:nRun
        if isstruct(models)
            allNames{iRun} = models(iRun).reg_names(:)';
        elseif iscell(models) && isstruct(models{iRun})
            allNames{iRun} = models{iRun}.reg_names(:)';
        elseif iscell(models) && iscell(models{iRun})
            allNames{iRun} = models{iRun}(:)';
        else
            error('models must be a struct array or cell array of models / name lists');
        end
        if size(stats{iRun}.beta, 1) ~= numel(allNames{iRun})
            error('Run %d: number of reg_names does not match rows of stats.beta', iRun);
        end
    end

    % --- Information regressors per run, and their union ---
    infoIdx   = cell(1, nRun);
    dataNames = {};
    for iRun = 1:nRun
        names = allNames{iRun};
        infoIdx{iRun} = find(~contains(names, nuisance_patterns));
        dataNames = union(dataNames, names(infoIdx{iRun}), 'stable');
    end
    nReg  = numel(dataNames);
    nChan = size(stats{1}.beta, 2);

    for iRun = 2:nRun
        if size(stats{iRun}.beta, 2) ~= nChan
            error('Run %d has a different number of channels', iRun);
        end
    end

    % --- Map each run's information columns into the union ---
    mapU = cell(1, nRun);
    for iRun = 1:nRun
        [~, mapU{iRun}] = ismember(allNames{iRun}(infoIdx{iRun}), dataNames);
        if numel(unique(mapU{iRun})) < numel(mapU{iRun})
            error('Run %d: duplicate information regressor names', iRun);
        end
    end

    % --- Inverse-covariance weighted combination, per channel ---
    beta_fe = zeros(nReg, nChan);
    covb_fe = nan(nReg, nReg, nChan);
    nEff = zeros(nReg, nChan);

    for iChan = 1:nChan
        Info = zeros(nReg);
        Wb   = zeros(nReg, 1);

        for iRun = 1:nRun
            ii = infoIdx{iRun};
            iu = mapU{iRun};

            b = stats{iRun}.beta(ii, iChan);
            C = stats{iRun}.covb(ii, ii, iChan);

            if any(~isfinite(b)) || any(~isfinite(C(:)))
                continue;   % skip runs with failed channel fits
            end

            P = pinv(C);
            Info(iu, iu) = Info(iu, iu) + P;
            Wb(iu)       = Wb(iu) + P * b;
            nEff(iu, iChan) = nEff(iu, iChan) + models(iRun).n_event(ii)';
        end

        ok = diag(Info) > 0;
        if any(ok)
            C_fe = pinv(Info(ok, ok));
            beta_fe(ok, iChan)     = C_fe * Wb(ok);
            covb_fe(ok, ok, iChan) = C_fe;
        end
        beta_fe(~ok, iChan) = NaN;
    end

    % --- Degrees of freedom (sum across runs) ---
    dfe = sum(cellfun(@(s) s.dfe, stats));

    % --- Standard errors, t and p ---
    se = nan(nReg, nChan);
    for iChan = 1:nChan
        se(:, iChan) = sqrt(diag(covb_fe(:, :, iChan)));
    end
    tstat = beta_fe ./ se;
    pval  = 2 * tcdf(-abs(tstat), dfe);

    % --- Output ---
    stats_fe          = struct();
    stats_fe.beta     = beta_fe;
    stats_fe.covb     = covb_fe;
    stats_fe.dfe      = dfe;
    stats_fe.tstat    = tstat;
    stats_fe.pval     = pval;
    stats_fe.reg_names = dataNames;
    stats_fe.nEff = nEff;

end