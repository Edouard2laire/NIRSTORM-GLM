function model = nst_glm_concatenate(modelA, modelB)

    if isempty(modelA) || isempty(modelB) 
        if isempty(modelA)
            model = modelB;
        elseif isempty(modelB) 
            model = modelA;
        else
            error('Both model cannot be empty')
        end
        return;
    end
    
    modelA.reg_names;
    modelB.reg_names;
    
    dt = modelA.time(2) - modelA.time(1);
    new_start = modelA.time(end) + dt;

    new_time =  [modelA.time , ...
                 modelB.time - modelB.time(1) + new_start];

    isNuisanceA = contains(modelA.reg_names, {'Constant', 'Poly'});
    isNuisanceB = contains(modelB.reg_names, {'Constant', 'Poly'});

    X_nuisance = blkdiag(modelA.X(:, isNuisanceA) , modelB.X(:, isNuisanceB));
    X_data = vertcat(modelA.X(:, ~isNuisanceA) , modelB.X(:, ~isNuisanceB));


    model       = nst_glm_initialize_model(new_time);
    model.X     = [X_data, X_nuisance];
    
    model.reg_names = [modelA.reg_names(~isNuisanceA), modelA.reg_names(isNuisanceA), modelB.reg_names(isNuisanceB)];
    model.hrf = modelA.hrf;
    model.accept_filter = [modelA.accept_filter(~isNuisanceA), modelA.accept_filter(isNuisanceA), modelB.accept_filter(isNuisanceB)];

end