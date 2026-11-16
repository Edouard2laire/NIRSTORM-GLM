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

    


end