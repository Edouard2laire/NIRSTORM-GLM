function T =  display_summary(fit)

    A = fit.summary();
    
    vals = double(A.to_numpy());
    rows = cellfun(@(x) string(x), cell(A.index.tolist()));
    cols = cellfun(@(x) string(x), cell(A.columns.tolist()));
    
    T = array2table(vals, ...
        'RowNames', rows, ...
        'VariableNames', matlab.lang.makeValidName(cols));
    
    if nargout < 1
        disp(T)
    end

end