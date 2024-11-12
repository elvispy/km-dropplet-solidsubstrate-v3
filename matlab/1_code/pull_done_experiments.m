% PULL_DONE_EXPERIMENTS Pulls data from completed simulation files into a table.
%
% Syntax:
%   A = pull_done_experiments(simulations)
%
% Description:
%   This function searches through subdirectories of the parent directory
%   for `.mat` files that contain simulation results. It filters out files
%   containing 'error' in their names, loads specific fields from each file,
%   and adds them as rows in an output table.
%
% Input:
%   simulations - (table) A table with column names corresponding to
%                 expected fields in the simulation `.mat` files.
%
% Output:
%   A - (table) A table of loaded simulations with rows corresponding to
%       each successfully loaded simulation file.
%
% Example:
%   A = pull_done_experiments(simulations);
%
% Notes:
%   The function assumes `.mat` files include two structures:
%   - `default_numerical`: Contains numerical parameters of the simulation.
%   - `default_physical`: Contains physical parameters of the simulation.
%
%   If a required field (`nu`) is missing from `default_physical`, it is
%   added with a value of 0. The modified file is saved.

function A = pull_done_experiments(simulations)
    % Initialize the output table A with variable names matching input table
    fldnames = simulations.Properties.VariableNames;
    A = array2table(double.empty(0, length(fldnames)), 'VariableNames', fldnames);
    
    % Locate all .mat files in parent subdirectories
    res = dir(fullfile(fileparts(fileparts(pwd)), '**/*.mat'));
    
    % Process each file, skipping those with 'error' in their name
    for jj = 1:length(res)
        if contains(lower(res(jj).name), 'error') 
            continue; 
        end
        A = addSimulation(A, res(jj)); 
    end
end

% ADDSIMULATION Adds a simulation file's data to the table if compatible.
%
% Syntax:
%   A = addSimulation(A, stru)
%
% Description:
%   Attempts to load a simulation file and merge its fields with table A.
%   If a required field is missing, it is initialized to a default value.
%
% Input:
%   A    - (table) The table to which new simulation data is added.
%   stru - (struct) Information about the `.mat` file, such as folder and name.
%
% Output:
%   A - (table) Updated table with added simulation data.
%
% Errors:
%   The function handles specific errors with custom messages, including:
%   - Size mismatch when the file's data structure does not match A.
%   - Undefined functions or variables if data is missing.
%
% Example:
%   A = addSimulation(A, res(jj));
%
function A = addSimulation(A, stru)
    try
        % Load simulation data
        load(fullfile(stru.folder, stru.name), 'default_numerical', 'default_physical');

        % Ensure 'nu' field is present in default_physical
        if any(strcmp(A.Properties.VariableNames, "nu")) && ...
                ~isfield(default_physical, 'nu')
            default_physical.nu = 0;
            save(fullfile(stru.folder, stru.name), "default_physical", "-append");
        end
        
        % Merge loaded structures
        tempStruct = cell2struct([struct2cell(default_numerical); struct2cell(default_physical)], ...
            [fieldnames(default_numerical); fieldnames(default_physical)]);
        
        % Keep only relevant fields in tempStruct
        tempStruct = rmfield(tempStruct, setdiff(fieldnames(tempStruct), A.Properties.VariableNames));
        
        % Append the structured data as a new row in A
        A = [A; struct2table(tempStruct)];
    catch me
        % Custom error handling for known issues
        switch me.identifier
            case 'MATLAB:table:vertcat:SizeMismatch'
                fprintf("Could not load simulation on %s, %s \n", stru.folder, stru.name);
            case 'MATLAB:UndefinedFunction'
                fprintf("Could not load simulation on %s, %s, undefined function or variable \n", stru.folder, stru.name);
            case 'MATLAB:load:unableToReadMatFile'
                warning("Possible corrupt file: %s/%s. It's going to be deleted \n",  stru.folder, stru.name);
                delete(fullfile(stru.folder, stru.name));
             otherwise
                rethrow(me)
        end
    end
end
