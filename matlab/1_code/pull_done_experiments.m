function A = pull_done_experiments(simulations)
    fldnames = simulations.Properties.VariableNames;
    A = array2table(double.empty(0, length(fldnames)), 'VariableNames', fldnames);
    
    res = dir("../2_output/*/*.mat");
    
    for jj = 1:length(res)
        if contains(res(jj).name, 'error') 
            continue; 
        end
        A = addSimulation(A, res(jj)); 
    end
    
   
end

function A = addSimulation(A, stru)
    load(fullfile(stru.folder, stru.name), 'default_numerical', 'default_physical');
    
    % Load structs and merge them
    tempStruct = cell2struct([struct2cell(default_numerical); struct2cell(default_physical)], ...
        [fieldnames(default_numerical); fieldnames(default_physical)]);
    
    
    try
        % Remove field that are not present 
        tempStruct = rmfield(tempStruct, setdiff(fieldnames(tempStruct), A.Properties.VariableNames));
        % Merge fields into simulared 
        A = [A;struct2table(tempStruct)];
    catch me
        switch me.identifier
            case 'MATLAB:table:vertcat:SizeMismatch'
                fprintf("Could not load simulation on %s, %s", stru.folder, stru.name);
            otherwise
                rethrow(me)
        end
    end
end