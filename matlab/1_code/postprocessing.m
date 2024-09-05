% Script to calculate relevant variables for all the simulations that
%have been run

% OS-independent parent folder
root_folder = fileparts(fileparts(mfilename('fullpath')));
% Setting diary path
diary(fullfile(root_folder, '0_data', 'manual', 'Logger', 'sweeper_postprocessing_logger.txt'));

disp("-------");
fprintf("%s \n %s", datestr(datetime()), mfilename('fullpath'));

% Add functions to calculate maximum width
safe_folder = fullfile(fileparts(mfilename('fullpath')), "simulation_code");
addpath(safe_folder, '-begin');

% FInding all .mat files that correspond to simulations
files_folder = dir(fullfile(root_folder, "2_output", "**/*.mat"));

% Create table to fill values
varNames=["fileName", "initialVelocity", "Westar", "Ohnesorge", "dropLiquid", "nb_harmonics", ...
    "maxWidth", "contactTime", "coefRestitution", "minHeight", "maxContactRadius"];
varTypes=["string", "double", "double", "double", "string", "double", "double", "double", "double", ...
    "double", "double"];
sz = [length(files_folder) length(varNames)];
data = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);


for ii = 1:length(files_folder)
    try
        lastwarn('', ''); clear recorded_conditions recorded_times default_physical length_unit theta_vector
        load(fullfile(files_folder(ii).folder, files_folder(ii).name), ...
            "recorded_conditions", "recorded_times", "default_physical", ...
            "length_unit", "PROBLEM_CONSTANTS");
        if contains(lastwarn, 'not found'); error(lastwarn);  end
        
        theta_vector = PROBLEM_CONSTANTS.theta_vector;
        Westar = default_physical.rhoS * default_physical.initial_velocity^2 * ...
            default_physical.undisturbed_radius / default_physical.sigmaS;
        Oh = default_physical.nu * sqrt(default_physical.rhoS / (default_physical.sigmaS ...
            * default_physical.undisturbed_radius));

        max_width = -inf;
        contact_time = nan; touch_time = nan; liftoff_time = nan;
        coef_restitution = nan; Vin = nan; Vout = nan;
        min_height = inf;
        max_contact_radius = -inf;
        for jj = 1:(size(recorded_conditions, 1)-1)
            adim_deformations = recorded_conditions{jj}.deformation_amplitudes/length_unit;
            adim_CM = recorded_conditions{jj}.center_of_mass/length_unit;
            drop_radius = zeta_generator(adim_deformations);
            drop_radius = @(theta) 1 + drop_radius(theta);

            % Max width calculation
            current_width = maximum_contact_radius(adim_deformations);
            if current_width > max_width; max_width = current_width; end

            % contact_time and coef_res
            if isnan(touch_time)
                if recorded_conditions{jj}.contact_points == 0 && ...
                        recorded_conditions{jj+1}.contact_points > 0
                    touch_time = recorded_times(jj);
                    Vin = recorded_conditions{jj}.center_of_mass_velocity;
                end
            end
            if isnan(liftoff_time)
                if recorded_conditions{jj}.contact_points > 0  && ...
                        recorded_conditions{jj+1}.contact_points == 0
                    liftoff_time = recorded_times(jj);
                    Vout = recorded_conditions{jj}.center_of_mass_velocity;
                end
            end
            if ~isnan(touch_time) && ~isnan(liftoff_time)
                contact_time = liftoff_time - touch_time;
                coef_restitution = abs(Vout/Vin);
            end

            % Min_height
            current_height = drop_radius(0) + adim_CM;
            if current_height < min_height; min_height = current_height; end

            % max_contact_radius calculation
            current_contact_points = recorded_conditions{jj}.contact_points;
            if current_contact_points == 0
                current_contact_radius = 0;
            else
                current_contact_radius = sin(theta_vector(current_contact_points)) ...
                    * drop_radius(theta_vector(current_contact_points));
            end
            if current_contact_radius > max_contact_radius; max_contact_radius = current_contact_radius; end
        end
        % Add value to tables
        [~, dropLiquid, lol] = fileparts(files_folder(ii).folder); dropLiquid = strcat(dropLiquid, lol);
        data(ii, :) = {files_folder(ii).name, Vin, Westar, Oh, dropLiquid, PROBLEM_CONSTANTS.nb_harmonics, ...
            max_width, contact_time, coef_restitution, min_height, max_contact_radius};
    catch me
        fprintf("%s %s \n", files_folder(ii).folder, files_folder(ii).name)
        disp(me);
    end
end
diary off
