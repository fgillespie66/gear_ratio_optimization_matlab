function pb = setParameterValues(pb)


fn = fieldnames(pb.p);
for j = 1:length(fn)
    try
        pb.opti.set_value(pb.p.(fn{j}).sym,pb.p.(fn{j}).val) 
    catch ME
        if strcmp(ME.identifier,'MATLAB:nonExistentField')
            msg = [ME.message ' for field name ', fn{j}];
            error(msg)
        else
            disp(['problem in setParameterValuesHybrid when setting value for field ', fn{j}]);
            rethrow(ME)
        end
    end
end

