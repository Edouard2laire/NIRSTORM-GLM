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


    model       = nst_glm_initialize_model(new_time);
    model.X     = blkdiag(modelA.X, modelB.X);
    
    model.reg_names = [modelA.reg_names, modelB.reg_names];
    model.hrf = modelA.hrf;
    model.accept_filter = [modelA.accept_filter, modelB.accept_filter];

end