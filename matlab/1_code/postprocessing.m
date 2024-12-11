% Script to calculate relevant variables for all the simulations that
%have been run

% OS-independent parent folder
root_folder = fileparts(fileparts(mfilename('fullpath')));
% Setting diary path
diary(fullfile(root_folder, '0_data', 'manual', 'Logger', 'sweeper_postprocessing_logger.txt'));

disp("-------");
fprintf("%s \n %s \n", datestr(datetime()), mfilename('fullpath'));

% Add functions to calculate maximum width
safe_folder = fullfile(fileparts(mfilename('fullpath')), "simulation_code");
addpath(safe_folder, '-begin');

% FInding all .mat files that correspond to simulations
files_folder = dir(fullfile(root_folder, "2_output", "**/*.mat"));


% Create table to fill values (adimensional unless stated in the title)
varNames=["file_name", "initial_velocity_cgs", "weber", "bond", "ohnesorge", "parent_folder", "number_of_harmonics", ...
    "max_width", "contact_time_ms", "coef_restitution", "north_pole_min_height", "north_pole_exp_min_height", ...
    "max_contact_radius", "spread_time_ms", "spread_time_width_ms", ...
    "contact_time_exp_ms", "coef_rest_exp", "max_contact_radius_exp", "spread_time_exp_ms"];
varTypes=["string", "double", "double", "double", "double", "string", "double", "double", "double", "double", ...
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
pixel = 5e-4; %Threshold for experimental contact
fnames = data{:, 1};
parfor ii = 1:length(files_folder)
    try
        if ismember(files_folder(ii).name, fnames) || ...
                contains(files_folder(ii).name, "postprocessing") || ...
                contains(lower(files_folder(ii).name), "error"); continue; 
        end
        %if ~isnan(data{ii, "coef_rest_exp"}); continue; end
        lastwarn('', ''); %clear recorded_conditions recorded_times default_physical length_unit theta_vector
        val = load(fullfile(files_folder(ii).folder, files_folder(ii).name), ...
            "recorded_conditions", "recorded_times", "default_physical", ...
            "length_unit", "PROBLEM_CONSTANTS");
        recorded_conditions = val.recorded_conditions;
        recorded_times = val.recorded_times;
        default_physical = val.default_physical;
        length_unit = val.length_unit;
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
        contact_time = nan; touch_time = nan; liftoff_time = nan;
        contact_time_exp = nan; liftoff_time_exp = nan;
        coef_restitution = nan; Vin = nan; Vout = nan;
        coef_restitution_exp = nan; Vout_exp = nan;
        north_pole_min_height = inf; north_pole_exp_min_height = inf;
        max_contact_radius = -inf; spread_time = nan; spread_time_width = nan;
        max_contact_radius_exp = -inf; spread_time_exp = nan; 
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
                    CM_in = recorded_conditions{jj}.center_of_mass;
                    Ein = 1/2 * Vin^2;
                end
            end
            
            if isnan(liftoff_time)
                if recorded_conditions{jj}.contact_points > 0  && ...
                        recorded_conditions{jj+1}.contact_points == 0
                    liftoff_time = recorded_times(jj);
                    Vout = recorded_conditions{jj}.center_of_mass_velocity;
                    Eout = 1/2*Vout^2 + (recorded_conditions{jj}.center_of_mass ...
                        - CM_in)*g;
                end
            end
            if ~isnan(touch_time) && ~isnan(liftoff_time)
                contact_time = liftoff_time - touch_time;
                coef_restitution = sqrt(abs(Eout/Ein));
            end

            % Min_height (north pole)
            current_height = drop_height(0);
            next_height = drop_height(pi/100);
            if current_height < north_pole_min_height; north_pole_min_height = current_height; end
            if next_height < current_height %&& current_height < north_pole_exp_min_height
                %north_pole_exp_min_height = current_height;
                if current_height < north_pole_exp_min_height; north_pole_exp_min_height = current_height; end
            else
                currang = pi/200; dth = pi/200;
                next_next_height = drop_height(currang + dth);
                while next_next_height >= next_height && currang < pi/3
                    currang = currang + dth;
                    next_height = next_next_height;
                    next_next_height = drop_height(currang + dth);
                end
                if next_height < north_pole_exp_min_height; north_pole_exp_min_height = next_height; end
            end

            % Experimental contact_radius
            if isnan(liftoff_time_exp)
                % If contact ended numerically but simulation ended, record lift off time anyways
                if drop_height(pi) > pixel/length_unit %||((size(recorded_conditions, 1)-1 == jj && ~isnan(liftoff_time)))
                    liftoff_time_exp = recorded_times(jj);
                    Vout_exp = recorded_conditions{jj}.center_of_mass_velocity;
                    Eout_exp = 1/2*Vout_exp^2 + (recorded_conditions{jj}.center_of_mass ...
                        - CM_in)*g;
                    % if (size(recorded_conditions, 1)-1 == jj && ~isnan(liftoff_time)) 
                    %     warning("Simulation ended abruptly but numerical contact has ended. We recorded the experimental end of contact anyways for simul \n %s", files_folder(ii).name);
                    % end
                    
                end
            end
                
            % Max width calculation & spread time of width
            current_width = maximum_contact_radius(adim_deformations);
            if current_width > max_width
                max_width = current_width; 
                spread_time_width = recorded_times(jj) - touch_time;
            end
            % We will record experimental touch time 
            if ~isnan(touch_time) && ~isnan(liftoff_time_exp) 
                contact_time_exp = liftoff_time_exp - touch_time;
                coef_restitution_exp = sqrt(abs(Eout_exp/Ein));
            end
            


            % max_contact_radius calculation
            current_contact_points = recorded_conditions{jj}.contact_points;
            if current_contact_points == 0
                current_contact_radius = 0;
            else
                current_contact_radius = sin(theta_vector(current_contact_points)) ...
                    * drop_radius(theta_vector(current_contact_points));
            end
            if current_contact_radius > max_contact_radius
                max_contact_radius = current_contact_radius;
                spread_time = recorded_times(jj) - touch_time;
            end

            % EXPERIMENTAL max_contact_radius calculation
            current_contact_points_exp = recorded_conditions{jj}.contact_points;
            if current_contact_points_exp == 0
                current_contact_radius_exp = 0;
            else
                
                kk = 1;
                theta2 = theta_vector(current_contact_points_exp+kk);
                while drop_height(theta2) < pixel/length_unit
                    kk = kk + 1;
                    theta2 = theta_vector(current_contact_points_exp + kk);
                end
                theta1 = theta_vector(current_contact_points_exp + kk - 1);
                while abs(theta1-theta2) > pi/1e+3
                    theta_middle = (theta1+theta2)/2;
                    if drop_height(theta_middle) < pixel/length_unit
                        theta1 = theta_middle;
                    else
                        theta2 = theta_middle;
                    end
                end
                current_contact_radius_exp = sin(theta1) ...
                    * drop_radius(theta1);
            end
            if current_contact_radius_exp > max_contact_radius_exp
                max_contact_radius_exp = current_contact_radius_exp;
                spread_time_exp = recorded_times(jj) - touch_time;
            end
        end
        % Add value to tables
        [~, dropLiquid, lol] = fileparts(files_folder(ii).folder); dropLiquid = strcat(dropLiquid, lol);
        data(ii, :) = {files_folder(ii).name, Vin, Westar, Bo, Oh, dropLiquid, PROBLEM_CONSTANTS.nb_harmonics, ...
            max_width, contact_time, coef_restitution, north_pole_min_height, ...
            north_pole_exp_min_height, max_contact_radius, spread_time, spread_time_width, ...
            contact_time_exp, coef_restitution_exp, max_contact_radius_exp, spread_time_exp};

    catch me
        if contains(files_folder(ii).name, "error"); continue; end
        fprintf("%s %s \n", files_folder(ii).folder, files_folder(ii).name)
        disp(me);
        switch me.identifier
            case 'MATLAB:load:unableToReadMatFile'
                delete(fullfile(iles_folder(ii).folder, files_folder(ii).name));
        end
    end
end
delete(gcp("nocreate")); % Deleting current parallel workers
% Filtering
data = rmmissing(data, 'DataVariables','file_name');

writetable(data, fullfile(root_folder, "2_output", "postprocessing.csv"));
s = fullfile(root_folder, "2_output", "postprocessing.mat");
if ~exist(s, "file")
    save(s, "data");
else
    warning('Did not save the postprocessing file. There is one in the directory already');
end
diary off
