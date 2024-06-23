% This sript will try to sweep simulations according to two
% rules:
% 1) Parameters set in this sweep
% 2) All empty folders will be swept. 

% STEP 1: Define which simulations are to be run. 

ddate = datestr(datetime(), 30); % 30 = ISO 8601
diary(sprintf("../2_output/Logger/%s_%s.txt", mfilename, ddate));
disp("------------");
fprintf("%s \n %s\n", string(datetime("now")), mfilename('fullpath'));


%% Setting simulation parameters
%#ok<*NOPTS>
vars = struct(...    %D = 50  %Quant = 100
    "RhoS", 1, ... % must multiply by x1000
    "SigmaS", 72.20, ... % must multiply by x100
    "R", 0.035, ... % linspace(0.02, 0.05, 5)'; % must multiply by x10 %Ang = 180
    "U", -linspace(10, 50, 5)', ... %inspace(59, 39, 6)';
    "modes", [10, 20, 40]', ...
    "version", [1, 2, 3]');%tol = 5e-5

% We check how many outputs we want
numOutputs = length(fieldnames(vars));


idxs = cell(1, numOutputs);
lengths = arrayfun(@(x) length(x{1}), struct2cell(vars), 'UniformOutput', false);
idxinputs = arrayfun(@(x) 1:x{1}, lengths, 'UniformOutput', false);
[idxs{:}] = ndgrid(idxinputs{:});

fnames = fieldnames(vars);
cartesian_product = cell2mat( ...
    arrayfun(@(varidx) vars.(fnames{varidx})(idxs{varidx}, :), 1:numOutputs, ...
    'UniformOutput',false));


% Turn simulations into table
if isempty(cartesian_product) == true; cartesian_product = double.empty(0, length(fieldnames(vars))); end
simulations_cgs = array2table(cartesian_product, "VariableNames", fnames);
% Now you can manually add any simulations that you would like to run, such
% as:
%  simulations_cgs = [simulations_cgs; ...
%      {50, 100, 1, 72.20, 0, 9.78E-3, 1, 72.20, 0.001, 180, 1}];
nrow = size(simulations_cgs, 1);
simulations_cgs.folder = repmat("../2_output/", nrow, 1);



% List of files needed to run the simulation and that do not need to change
% from simul to simul.

force_sweep = false;

% STEP 3: Actually run the simulations. 

% Initial folder to go back to
root = pwd;
% A folder which MUST have all the dependencies needed (and all its
% parents, too. 
safe_folder = fullfile(root, "simulation_code");
addpath(safe_folder, '-begin');
%safe_folder = fullfile(root, "D50Quant100", "rho1000sigma7220nu98muair0", ...
%    "RhoS1000SigmaS7220", "R0350mm", "ImpDefCornerAng180U38");

final_folders = simulations_cgs.folder;

% Some common features of all the simulations to be run


%% Starting simulation
parfor ii = 1:height(simulations_cgs)
    %Check if etaOri exists (the center of the bath)
    cd(final_folders(ii));

    if force_sweep == true || isempty(dir("recorded*.mat")) == true
        
        numerical_parameters = struct("harmonics_qtt", simulations_cgs.modes(ii),...
            "simulation_time", inf, "version", simulations_cgs.version(ii));
        options = struct('version', simulations_cgs.version(ii), ...
            'folder', sprintf("Version v%d (%g)", simulations_cgs.version(ii), ...
        simulations_cgs.RhoS(ii) * simulations_cgs.SigmaS(ii)));

        physical_parameters = struct("undisturbed_radius", simulations_cgs.R(ii), "initial_height", nan, ...
            "initial_velocity", simulations_cgs.U(ii), "initial_amplitudes", zeros(1, simulations_cgs.modes(ii)), ...
            "pressure_amplitudes", zeros(1, simulations_cgs.modes(ii)+1), "initial_contact_points", 0, ...
            "rhoS", simulations_cgs.RhoS(ii), "sigmaS", simulations_cgs.SigmaS(ii));
        
        %solve_motion_v2(physical_parameters, numerical_parameters);

        try
            fprintf("Starting simulation with velocity %g, modes %d, version v%d", ...
                simulations_cgs.U(ii), simulations_cgs.modes(ii), simulations_cgs.version(ii));
            solve_motion_v2(physical_parameters, numerical_parameters, options);
        catch ME
            cd(final_folders(ii))
            fprintf("Couldn't run simulation with the following parameters: \n Velocity: %g \n Modes: %g \n", ...
                simulations_cgs.U(ii), simulations_cgs.modes(ii)); 
            a = datetime('now'); a.Format = 'yyyyMMddmmss';
            parsave(sprintf("error_logU0=%g-%s.mat", simulations_cgs.U(ii), a), ME);
        end
    else
        fprintf("Not running simulation with the following parameters (already done): \n Velocity: %g \n Modes: %g \n", ...
                simulations_cgs.U(ii), simulations_cgs.modes(ii)); 
    end
    

end

cd(root);
delete(gcp("nocreate")); % Deleting current parallel workers

% Load Python3 in MACOS based on https://www.mathworks.com/matlabcentral/answers/359408-trouble-with-a-command-in-matlab-s-system
if ~ispc && system('python3 --version') ~= 0; setenv('PATH', [getenv('PATH') ':/usr/local/bin/']); end

false && system('python3 sending_email.py'); % Sending email to notify that's finished
diary off % turning logger off


