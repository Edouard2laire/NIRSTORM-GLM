function model = nst_glm_concatenate(modelA, modelB)

    % --- Handle empty inputs ---
    if isempty(modelA) && isempty(modelB)
        model = modelB;   % both empty: return empty
        return;
    elseif isempty(modelA)
        model = modelB;
        return;
    elseif isempty(modelB)
        model = modelA;
        return;
    end

    % Force row cell arrays / row logical vectors for consistent indexing
    namesA = modelA.reg_names(:)';
    namesB = modelB.reg_names(:)';
    filtA  = modelA.accept_filter(:)';
    filtB  = modelB.accept_filter(:)';

    % --- New time vector ---
    dt        = modelA.time(2) - modelA.time(1);
    new_start = modelA.time(end) + dt;
    new_time  = [modelA.time(:)', ...
                 modelB.time(:)' - modelB.time(1) + new_start];

    % --- Split nuisance / information regressors ---
    isNuisanceA = contains(namesA, {'Constant', 'Poly'});
    isNuisanceB = contains(namesB, {'Constant', 'Poly'});

    dataNamesA = namesA(~isNuisanceA);
    dataNamesB = namesB(~isNuisanceB);

    if numel(unique(dataNamesA)) < numel(dataNamesA) || ...
       numel(unique(dataNamesB)) < numel(dataNamesB)
        error('Duplicate information regressor names found within a model');
    end

    XdataA = modelA.X(:, ~isNuisanceA);
    XdataB = modelB.X(:, ~isNuisanceB);
    filtDataA = filtA(~isNuisanceA);
    filtDataB = filtB(~isNuisanceB);

    % --- Union of information regressors (A order first, then new ones from B) ---
    data_names = union(dataNamesA, dataNamesB, 'stable');
    nReg = numel(data_names);
    nA   = size(modelA.X, 1);
    nB   = size(modelB.X, 1);

    X_data = zeros(nA + nB, nReg);
    data_filter = false(1, nReg);

    [inA, idxA] = ismember(data_names, dataNamesA);
    [inB, idxB] = ismember(data_names, dataNamesB);

    X_data(1:nA,     inA) = XdataA(:, idxA(inA));
    X_data(nA+1:end, inB) = XdataB(:, idxB(inB));

    % Accept filter: if present in both, require both to accept;
    % if present in only one, use that one's value
    fA = true(1, nReg);  fA(inA) = logical(filtDataA(idxA(inA)));
    fB = true(1, nReg);  fB(inB) = logical(filtDataB(idxB(inB)));
    data_filter = fA & fB;

    % --- Nuisance regressors: one set per file (block diagonal) ---
    X_nuisance = blkdiag(modelA.X(:, isNuisanceA), modelB.X(:, isNuisanceB));

    % --- Assemble ---
    model            = nst_glm_initialize_model(new_time);
    model.X          = [X_data, X_nuisance];
    model.reg_names  = [data_names, namesA(isNuisanceA), cellfun(@(x) file_unique(x, namesA(isNuisanceA)), namesB(isNuisanceB), 'UniformOutput', false)];
    model.hrf        = modelA.hrf;
    model.accept_filter = [data_filter, filtA(isNuisanceA), filtB(isNuisanceB)];
end