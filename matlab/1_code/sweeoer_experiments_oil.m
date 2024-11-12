% This sript will try to sweep simulations according to two
% rules:
% 1) Parameters set in this sweep
% 2) All empty folders will be swept. 

% STEP 1: Define which simulations are to be run. 

ddate = datestr(datetime(), 30); % 30 = ISO 8601
diary(sprintf("../2_output/Logger/%s_%s.txt", mfilename, ddate));
disp("------------");
fprintf("%s \n %s\n", string(datetime("now")), mfilename('fullpath'));
force_sweep = true;

%% Setting simulation parameters
%#ok<*NOPTS>
prefix = 'WeBoSweep';
sigma = 20.5; rho = 0.96; Ro = 0.0203;
We = logspace(-5, 0, 41); %10.^([-3, -2, -1, 0]);
Bo = 10.^([-inf, -3, -2, -1, 0]);
velocities = -sqrt(sigma/(rho*Ro) .* We); % [V, 2*V 3*V, 4*V, 5*V, 6*V, 7*V, 8*V, 9*V]);
g = Bo.* sigma ./(rho .* Ro^2);
vars = struct(...  
    "rhoS", rho, ... % Droplet density in cgs
    "sigmaS", sigma, ... % Surface tension in cgs
    "nu", 0.01 * [1, 20]', ... % Viscocity in cgs
    "g", g', ... % Gravity (cgs)
    "undisturbed_radius", Ro, ... % (cgs)
    "initial_velocity", velocities', ... %(cgs)
    "harmonics_qtt", [90]', ...
    "version", [3]')

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

% retrieve simulations that have already been done. 
if force_sweep == false
    simulations_cgs.done = ismember(simulations_cgs, pull_done_experiments(simulations_cgs)); 
else
    simulations_cgs.done = zeros(height(simulations_cgs), 1);
end

root = pwd;
parentDir = fileparts(fileparts(mfilename('fullpath')));

a=rowfun(@(a, b, c, d) fullfile(parentDir, '2_output', sprintf("Version v%d (rhoS=%.2e, sigmaS=%.2e, R=%.2e)", a, b, c, d)), ...
    simulations_cgs, 'InputVariables', {'version' 'rhoS', 'sigmaS', 'undisturbed_radius'}, 'OutputVariableName', 'folder');
simulations_cgs.folder = a.folder;

% STEP 3: Actually run the simulations. 


% A folder which MUST have all the dependencies needed (and all its
% parents, too. 
safe_folder = fullfile(root, "simulation_code");
addpath(safe_folder, '-begin');
%safe_folder = fullfile(root, "D50Quant100", "rho1000sigma7220nu98muair0", ...
%    "RhoS1000SigmaS7220", "R0350mm", "ImpDefCornerAng180U38");

%final_folders = simulations_cgs.folder;

% Extract columns from table as parfor does not allowusing them inside the
% loop (TOOD: convert the table to a struct array
final_folders = simulations_cgs.folder;
completed_simulations = simulations_cgs.done;
harmonics_qtt = simulations_cgs.harmonics_qtt;
version = simulations_cgs.version;
rhoS = simulations_cgs.rhoS;
sigmaS = simulations_cgs.sigmaS;
nu = simulations_cgs.nu; g = simulations_cgs.g;
initial_velocity = simulations_cgs.initial_velocity;
undisturbed_radius = simulations_cgs.undisturbed_radius;
%% Starting simulation
parfor ii = 1:height(simulations_cgs)
    %Check if the final folder exists
    if ~exist(final_folders(ii), 'dir')
        mkdir(final_folders(ii))
    else
        cd(final_folders(ii));
    end

    if completed_simulations(ii) == false
        
        numerical_parameters = struct("harmonics_qtt", harmonics_qtt(ii),...
            "simulation_time", inf, "version", version(ii));
        options = struct('version', version(ii), 'prefix', prefix);

        physical_parameters = struct("undisturbed_radius", undisturbed_radius(ii), ...
            "initial_height", nan, "initial_velocity", initial_velocity(ii), ...
            "initial_amplitudes", zeros(1, harmonics_qtt(ii)), ...
            "pressure_amplitudes", zeros(1, harmonics_qtt(ii)+1), "initial_contact_points", 0, ...
            "rhoS", rhoS(ii), "sigmaS", sigmaS(ii), "nu", nu(ii), "g", g(ii));
        
        %solve_motion_v2(physical_parameters, numerical_parameters);

        try
            fprintf("---------\n");
            fprintf("Starting simulation with velocity %g, modes %d, version v%d ... \n", ...
                initial_velocity(ii), harmonics_qtt(ii), version(ii));
            solve_motion_v2(physical_parameters, numerical_parameters, options);
            completed_simulations(ii) = true; % To attest that the simulation has been finished
            fprintf("Finished simulation with velocity %g, modes %d, version v%d ... \n", ...
                initial_velocity(ii), harmonics_qtt(ii), version(ii));
        catch ME
            cd(final_folders(ii))
            fprintf("---------\n");
            fprintf("Couldn't run simulation with the following parameters: \n Velocity: %g \n Modes: %g \n Version: %g \n", ...
                initial_velocity(ii), harmonics_qtt(ii), version(ii)); 
            a = datetime('now'); a.Format = 'yyyyMMddmmss';
            parsave(sprintf("error_logU0=%g-%s.mat", initial_velocity(ii), a), ME);
        end
    else
        fprintf("---------\n");
        fprintf("Not running simulation with the following parameters (already done): \n Velocity: %g \n Modes: %g \n Version: %g \n", ...
                initial_velocity(ii), harmonics_qtt(ii), version(ii)); 
    end
    

end

cd(root);
delete(gcp("nocreate")); % Deleting current parallel workers

% Load Python3 in MACOS based on https://www.mathworks.com/matlabcentral/answers/359408-trouble-with-a-command-in-matlab-s-system
if ~ispc && system('python3 --version') ~= 0; setenv('PATH', [getenv('PATH') ':/usr/local/bin/']); end

system('python3 sending_email.py'); % Sending email to notify that's finished
diary off % turning logger off


function parsave(fname, errormsg)
    save(fname, 'errormsg')
end
