% This sript will try to sweep simulations according to two
% rules:
% 1) Parameters set in this sweep
% 2) All empty folders will be swept. 

% STEP 1: Define which simulations are to be run. 

ddate = datestr(datetime(), 30); % 30 = ISO 8601
diary(sprintf("../2_output/Logger/%s_%s.txt", mfilename, ddate));
disp("------------");
fprintf("%s \n %s\n", string(datetime("now")), mfilename('fullpath'));
force_sweep = false;
postprocessing_bool = true;
optimize_for_bounce = true;


%% Setting simulation parameters
%#ok<*NOPTS>
prefix = 'energyPlot';
sigma = 20.5; rho = 0.96; Ro = 0.0203;
We_exp = [0.254, 0.250] + rand()*1e-4;
Bo = 0.0189;%10.^([-inf, -3, -2, -1, 0]);
velocities_exp = -sqrt(sigma/(rho*Ro) .* We_exp); % [V, 2*V 3*V, 4*V, 5*V, 6*V, 7*V, 8*V, 9*V]);
g = Bo.* sigma ./(rho .* Ro^2);
X =2*g*Ro*0.02;
velocities = -sqrt(velocities_exp.^2 + X);
vars = struct(...  
    "rhoS", rho, ... % Droplet density in cgs
    "sigmaS", sigma, ... % Surface tension in cgs
    "nu", 0.01 * [2]', ... % Viscocity in cgs
    "g", g', ... % Gravity (cgs)
    "undisturbed_radius", Ro, ... % (cgs)
    "initial_velocity", velocities', ... %(cgs)
    "harmonics_qtt", [90, 150]', ...
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


files_folder = dir(fullfile(pwd, "2_output", "**/criticalWeOh180*.mat"));
We = zeros(1, length(files_folder));
Oh = zeros(1, length(files_folder));
Bo = zeros(1, length(files_folder));
new_simuls = zeros(length(files_folder), length(fnames));

% Sort by age (numeric field)
for ii = 1:length(files_folder)


    %if ~isnan(data{ii, "coef_rest_exp"}); continue; end
    lastwarn('', ''); %clear recorded_conditions recorded_times default_physical length_unit theta_vector
    % Clear unused variables from previous loops. 
    val = load(fullfile(files_folder(ii).folder, files_folder(ii).name), ...
        "recorded_conditions", "recorded_times", "default_physical", ...
        "length_unit", "PROBLEM_CONSTANTS", "time_unit");
    recorded_conditions = val.recorded_conditions;
    recorded_times = val.recorded_times;
    default_physical = val.default_physical;
    length_unit = val.length_unit;
    time_unit = val.time_unit;
    PROBLEM_CONSTANTS = val.PROBLEM_CONSTANTS;


    new_simuls(ii, :) = [default_physical.rhoS, default_physical.sigmaS, default_physical.nu, ...
        default_physical.g, default_physical.undisturbed_radius, default_physical.initial_velocity, ...
        90, 3]; % LAST TWO: HARMONICS_QTT AND VERSION;
end


% Now you can manually add any simulations that you would like to run, such
% as:
%  simulations_cgs = [simulations_cgs; ...
%      {50, 100, 1, 72.20, 0, 9.78E-3, 1, 72.20, 0.001, 180, 1}];

simulations_cgs = [simulations_cgs; array2table(new_simuls, 'VariableNames', simulations_cgs.Properties.VariableNames)];

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
%    "RhoS1000Sigma_folders = simulations_cgs.folder;

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
%% Starting simulatioS7220", "R0350mm", "ImpDefCornerAng180U38");

%finaln
T = height(simulations_cgs);
for ii = 1:height(simulations_cgs)
    %Check if the final folder exists
    if ~exist(final_folders(ii), 'dir')
        mkdir(final_folders(ii))
    else
        cd(final_folders(ii));
    end

    if completed_simulations(ii) == false
        
        numerical_parameters = struct("harmonics_qtt", harmonics_qtt(ii),...
            "simulation_time", inf, "version", version(ii));
        options = struct('version', version(ii), 'prefix', prefix, 'live_plotting', false,...
            'optimize_for_bounce', optimize_for_bounce);

        physical_parameters = struct("undisturbed_radius", undisturbed_radius(ii), ...
            "initial_height", nan, "initial_velocity", initial_velocity(ii), ...
            "initial_amplitudes", zeros(1, harmonics_qtt(ii)), ...
            "pressure_amplitudes", zeros(1, harmonics_qtt(ii)+1), "initial_contact_points", 0, ...
            "rhoS", rhoS(ii), "sigmaS", sigmaS(ii), "nu", nu(ii), "g", g(ii));
        
        %solve_motion_v2(physical_parameters, numerical_parameters);

        try
            %fprintf("---------\n");
            fprintf("Starting simulation %d/%d, with velocity %g, modes %d, version v%d ... \n", ...
                ii, T, initial_velocity(ii), harmonics_qtt(ii), version(ii));
            solve_motion_v2(physical_parameters, numerical_parameters, options);
            completed_simulations(ii) = true; % To attest that the simulation has been finished
            fprintf("Finished simulation %d/%d, with velocity %g, modes %d, version v%d ... \n", ...
                ii, T, initial_velocity(ii), harmonics_qtt(ii), version(ii));
        catch ME
            cd(final_folders(ii))
            fprintf("---------\n");
            fprintf("Couldn't run simulation %d/%d, with the following parameters: \n Velocity: %g \n Modes: %g \n Version: %g \n", ...
                ii, T, initial_velocity(ii), harmonics_qtt(ii), version(ii)); 
            a = datetime('now'); a.Format = 'yyyyMMddmmss';
            parsave(sprintf("error_logU0=%g-%s.mat", initial_velocity(ii), a), ME);
        end
    else
        fprintf("---------\n");
        fprintf("Not running simulation %d/%d, with the following parameters (already done): \n Velocity: %g \n Modes: %g \n Version: %g \n", ...
                ii, T, initial_velocity(ii), harmonics_qtt(ii), version(ii)); 
    end
    

end

cd(root);
% Load Python3 in MACOS based on https://www.mathworks.com/matlabcentral/answers/359408-trouble-with-a-command-in-matlab-s-system
if ~ispc && system('python3 --version') ~= 0; setenv('PATH', [getenv('PATH') ':/usr/local/bin/']); end

if postprocessing_bool == true
    postprocessing; 
else
    system('python3 sending_email.py'); % Sending email to notify that's finished
end
delete(gcp("nocreate")); % Deleting current parallel workers



diary off % turning logger off


function parsave(fname, errormsg)
    save(fname, 'errormsg')
end
