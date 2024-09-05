function A = pull_done_experiments(simulations)
    fldnames = simulations.Properties.VariableNames;
    A = array2table(double.empty(0, length(fldnames)), 'VariableNames', fldnames);
    
    res = dir(fullfile(fileparts(fileparts(pwd)), '**/*.mat')); %dir("../2_output/*/*.mat");
    
    for jj = 1:length(res)
        if contains(res(jj).name, 'error') 
            continue; 
        end
        A = addSimulation(A, res(jj)); 
    end
    
   
end

function A = addSimulation(A, stru)
	try
        load(fullfile(stru.folder, stru.name), 'default_numerical', 'default_physical');

        %Viscocity is special: it was added afterwards: If no viscocity is
        %present, we set it to zero
        if any(strcmp(A.Properties.VariableNames, "nu")) && ...
                ~isfield(default_physical, 'nu')
            default_physical.nu = 0;
            save(fullfile(stru.folder, stru.name), "default_physical", "-append");
        end
        
         % Load structs and merge them
        tempStruct = cell2struct([struct2cell(default_numerical); struct2cell(default_physical)], ...
            [fieldnames(default_numerical); fieldnames(default_physical)]);
        
        % Remove field that are not present 
        tempStruct = rmfield(tempStruct, setdiff(fieldnames(tempStruct), A.Properties.VariableNames));
        % Merge fields into simulared 
        A = [A;struct2table(tempStruct)];
    catch me
        switch me.identifier
            case 'MATLAB:table:vertcat:SizeMismatch'
                fprintf("Could not load simulation on %s, %s \n", stru.folder, stru.name);
            case 'MATLAB:UndefinedFunction'
                fprintf("Could not load simulation on %s, %s, undefined function or variable \n", stru.folder, stru.name);
            otherwise
                rethrow(me)
        end % end switch
    end % end try
end