% Script to calculate relevant variables for all the simulations that
%have been run. This script will calcualte only relevant variables 
% With two tresholds: 0.00R and 0.02R. 

% OS-independent parent folder
root_folder = fileparts(fileparts(mfilename('fullpath')));
% Setting diary path
diary(fullfile(root_folder, '0_data', 'manual', 'Logger', 'sweeper_postprocessing_sensitivity_logger.txt'));

disp("-------");
fprintf("%s \n %s \n", datestr(datetime()), mfilename('fullpath'));

% Add functions to calculate maximum width
safe_folder = fullfile(fileparts(mfilename('fullpath')), "simulation_code");
addpath(safe_folder, '-begin');

% FInding all .mat files that correspond to simulations
files_folder = dir(fullfile(root_folder, "2_output", "**/*.mat"));


% Create table to fill values (adimensional unless stated in the title)
varNames=["file_name", "initial_velocity_cgs", "weber", "bond", "ohnesorge", "parent_folder", "number_of_harmonics", ...
    "contact_time_ms_00R", "contact_time_ms002R",  "coef_restitution00R", "coef_restitution002R", ...
    "max_contact_radius_cm_00R", "max_contact_radius_cm_002R", "min_height_exp", "spread_time_ms_00R", "spread_time_ms_02R", ...
    "spread_time_width_ms"];
varTypes=["string", "double", "double", "double", "double", "string", "double", "double", ...
    "double", "double", "double", "double", "double", "double", "double", "double", "double"];

sz = [length(files_folder) length(varNames)];

if isfile(fullfile(root_folder, "2_output", "postprocessing.mat"))
    load(fullfile(root_folder,"2_output", "postprocessing.mat"), "data");
else
    data = table();
end

if size(data, 2) ~= length(varNames) || sum(data.Properties.VariableNames ~= varNames) > 2
    data = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);
end

data = data(~isnan(data.initial_velocity_cgs), :);
prefix = "higherResolution";
%data = data(contains(data.file_name, prefix), :);

%pixel = 5e-4; %Threshold for experimental contact
fnames = data{:, 1};

parfor ii = 1:length(files_folder)
    try
        if ismember(files_folder(ii).name, fnames) || ...
                ~contains(files_folder(ii).name, prefix) || ...
                contains(lower(files_folder(ii).name), "error"); continue; 
        end
        %if ~contains(files_folder(ii).name, "HighWeSweep20241110T150017"); continue; end
        
        %if ~isnan(data{ii, "coef_rest_exp"}); continue; end
        lastwarn('', ''); %clear recorded_conditions recorded_times default_physical length_unit theta_vector
        val = load(fullfile(files_folder(ii).folder, files_folder(ii).name), ...
            "recorded_conditions", "recorded_times", "default_physical", ...
            "length_unit", "PROBLEM_CONSTANTS", "time_unit");
        recorded_conditions = val.recorded_conditions;
        recorded_times = val.recorded_times;
        default_physical = val.default_physical;
        length_unit = val.length_unit;
        time_unit = val.time_unit;
        PROBLEM_CONSTANTS = val.PROBLEM_CONSTANTS;

        if contains(lastwarn, 'not found'); error(lastwarn);  end
        
        recorded_times = recorded_times * 1e+3; % To miliseconds
        sigmaS = default_physical.sigmaS; rhoS = default_physical.rhoS; Ro = default_physical.undisturbed_radius;
        Westar = rhoS * default_physical.initial_velocity^2 * Ro/sigmaS;
        
        if isfield(default_physical, 'nu'); Oh = default_physical.nu / sqrt(sigmaS * Ro * rhoS); else; Oh = 0; end
        g = default_physical.g;
        theta_vector = val.PROBLEM_CONSTANTS.theta_vector;
        Westar = default_physical.rhoS * default_physical.initial_velocity^2 * ...
            default_physical.undisturbed_radius / default_physical.sigmaS;
        Oh = default_physical.nu * sqrt(default_physical.rhoS / (default_physical.sigmaS ...
            * default_physical.undisturbed_radius));
        Bo = rhoS * default_physical.g * Ro^2 / sigmaS;
        
        pixel = 0.02 * length_unit;
        max_width = -inf;
        contact_time00R = nan; touch_time = nan; liftoff_time = nan;
        V0 = abs(default_physical.initial_velocity);
        t0 = (-V0 + sqrt(V0^2 - 2*g*0.02*Ro))/g; % Calculating experimental start of contact to substract
        if isnan(t0); t0 = -0.02*Ro/V0; end
        V = @(t)  -V0 - g*t; % COM Velocity
        X = @(t) -t * V0 - g/2*t.^2 +Ro; % COM position
        contact_time02R = nan; liftoff_time_exp = nan;
        coef_restitution00R = nan; Vin = nan; Vout = nan;
        coef_restitution02R = nan; Vout_exp = nan;
        %north_pole_min_height = inf; north_pole_exp_min_height = inf;
        max_contact_radius00R = -inf; spread_time_ms_00R = nan; spread_time_ms_02R = nan;
        max_contact_radius02R = -inf; spread_time_width_ms = nan; 
        north_pole_min_height = inf;
        
        velocity_unit = length_unit/time_unit;
        %energy_constant = sigmaS/(rhoS * Ro * velocity_unit.^2);
        %idxs = 1:PROBLEM_CONSTANTS.nb_harmonics;
        %Xl = (2*pi./(idxs .* (2 * idxs + 1))); Yl = (2*pi * (idxs.^2 + idxs - 2)./(2*idxs+1));
        %deformation_modes_energies = nan;
        %A = (velocity_unit.^2/default_physical.initial_velocity.^2)/(2*pi/3);

        for jj = 1:(size(recorded_conditions, 1)-1)
            adim_deformations = recorded_conditions{jj}.deformation_amplitudes/length_unit;
            adim_CM = recorded_conditions{jj}.center_of_mass/length_unit;
            drop_radius = zeta_generator(adim_deformations);
            drop_radius = @(theta) 1 + drop_radius(theta);
            drop_height = @(theta) cos(theta) .* drop_radius(theta) + adim_CM;

            
            % contact_time and coef_res
            if isnan(touch_time)
                if recorded_conditions{jj}.contact_points == 0 && ...
                        recorded_conditions{jj+1}.contact_points > 0
                    touch_time = recorded_times(jj);
                    Vin = recorded_conditions{jj}.center_of_mass_velocity;
                    CM_in    = recorded_conditions{jj}.center_of_mass; % 0.00R
                    CM_in02R = X(t0); % 0.02R
                    Ein    = 1/2 * Vin^2;  % 0.00R
                    Ein02R = 1/2 * V(t0)^2;% 0.02R
                end
            end
            
            if isnan(liftoff_time)
                if recorded_conditions{jj}.contact_points > 0  && ...
                        recorded_conditions{jj+1}.contact_points == 0
                    liftoff_time = recorded_times(jj);
                    Vout = recorded_conditions{jj}.center_of_mass_velocity; % 0.00R
                    Eout = 1/2*Vout^2 + (recorded_conditions{jj}.center_of_mass ...
                        - CM_in)*g; %0.00R
                end
            end
            if ~isnan(touch_time) && ~isnan(liftoff_time)
                contact_time00R = liftoff_time - touch_time;
                coef_restitution00R = sqrt(abs(Eout/Ein));
            end

            % Min_height (north pole)
            north_pole_exp_max_projection = drop_height(0); currang = pi/4; dtheta = pi/2; N = 9;
            while dtheta > pi/N^3
               vals = drop_height(linspace(currang-dtheta/2, currang+dtheta/2, N));
               [maxval, maxindex] = max(vals);
               if maxval > north_pole_exp_max_projection
                   north_pole_exp_max_projection = maxval;
               end
               currang = max(currang - dtheta/2 + (maxindex-1)/N * dtheta, dtheta/(2*N));
               dtheta = dtheta/N;
               if north_pole_exp_max_projection == drop_height(0)
                   break;
               end
            end
            
            north_pole_min_height = min(north_pole_min_height, north_pole_exp_max_projection);

            % Experimental contact_radius && coef restitution
            if isnan(liftoff_time_exp)
                % If contact ended numerically but simulation ended, record lift off time anyways
                if drop_height(pi) > pixel/length_unit %||((size(recorded_conditions, 1)-1 == jj && ~isnan(liftoff_time)))
                    liftoff_time_exp = recorded_times(jj);
                    Vout_exp = recorded_conditions{jj}.center_of_mass_velocity; % 0.02R
                    Eout_exp = 1/2*Vout_exp^2 + (recorded_conditions{jj}.center_of_mass ...
                        - CM_in02R)*g; % 0.02R

                end
            end
                
            % % Max width calculation & spread time of width
            current_width = maximum_contact_radius(adim_deformations);
            if current_width > max_width
                max_width = current_width; 
                spread_time_width_ms = recorded_times(jj) - touch_time;
            end
            % We will record experimental contact time
            if ~isnan(touch_time) && ~isnan(liftoff_time_exp) 
                contact_time02R = liftoff_time_exp - (touch_time + t0);
                coef_restitution02R = sqrt(abs(Eout_exp/Ein02R)); %  0.02R
            end
            
            % max_contact_radius calculation
            current_contact_points = recorded_conditions{jj}.contact_points;
            kk = 10;
            theta2 = theta_vector(current_contact_points+kk);
            while drop_height(theta2) < 1e-6
                kk = kk + 1;
                theta2 = theta_vector(current_contact_points + kk);
            end
            theta1 = theta_vector(current_contact_points + kk - 1);
            while abs(theta1-theta2) > pi/1e+4
                theta_middle = (theta1+theta2)/2;
                if drop_height(theta_middle) < 1e-6
                    theta1 = theta_middle;
                else
                    theta2 = theta_middle;
                end
            end
            current_contact_radius = sin(theta_middle) ...
                * drop_radius(theta_middle) * length_unit;
            
            if current_contact_radius > max_contact_radius00R
                max_contact_radius00R = current_contact_radius; %0.00R
                spread_time_ms_00R = recorded_times(jj) - touch_time;
            end

            % EXPERIMENTAL max_contact_radius calculation
            current_contact_points_exp = recorded_conditions{jj}.contact_points;
            kk = 1;
            theta2 = theta_vector(current_contact_points_exp+kk);
            while drop_height(theta2) < pixel/length_unit
                kk = kk + 1;
                theta2 = theta_vector(current_contact_points_exp + kk);
            end
            theta1 = theta_vector(current_contact_points_exp + kk - 1);
            while abs(theta1-theta2) > pi/1e+4
                theta_middle = (theta1+theta2)/2;
                if drop_height(theta_middle) < pixel/length_unit
                    theta1 = theta_middle;
                else
                    theta2 = theta_middle;
                end
            end
            current_contact_radius_exp = sin(theta1) ...
                * drop_radius(theta1) * length_unit;

            if current_contact_radius_exp > max_contact_radius02R
                max_contact_radius02R = current_contact_radius_exp; % 0.02R
                spread_time_ms_02R = recorded_times(jj) - touch_time;
            end
        end
        % Add value to tables
        [~, dropLiquid, lol] = fileparts(files_folder(ii).folder); dropLiquid = strcat(dropLiquid, lol);
        data(ii, :) = {files_folder(ii).name, Vin, Westar, Bo, Oh, dropLiquid, PROBLEM_CONSTANTS.nb_harmonics, ...
            contact_time00R, contact_time02R, coef_restitution00R, coef_restitution02R, ...
            max_contact_radius00R, max_contact_radius02R, north_pole_min_height, ...
            spread_time_ms_00R, spread_time_ms_02R, spread_time_width_ms};

    catch me
        if contains(files_folder(ii).name, "error"); continue; end
        fprintf("Exception at: %s %s \n", files_folder(ii).folder, files_folder(ii).name)
        disp(me);
        switch me.identifier
            case 'MATLAB:load:unableToReadMatFile'
                disp(" deletingfiles!");
                delete(fullfile(iles_folder(ii).folder, files_folder(ii).name));
        end
    end
end
delete(gcp("nocreate")); % Deleting current parallel workers
% Filtering
data = rmmissing(data, 'DataVariables','file_name');
%data_chase = data(data.number_of_harmonics == 90 & contains(data.parent_folder, 'v3') & abs(data.ohnesorge - 0.0303767)./data.ohnesorge <= 0.01 & data.bond <= 0.0191, :);
%writetable(data_chase, fullfile(root_folder, "2_output", "data_chase.csv"));
s = fullfile(root_folder, "2_output", "postprocessing_sensitivity.mat");
warning ('on','all');
if ~exist(s, "file")
    save(s, "data");
else
    warning('Did not save the postprocessing file. There is one in the directory already');
end
diary off
