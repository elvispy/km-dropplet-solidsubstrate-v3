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
if exist('prefix', 'var') ~= 1; prefix = ""; end % Only look for simulations starting with this prefix
files_folder = dir(fullfile(root_folder, "2_output", sprintf("**/%s*.mat", prefix)));

debug = false;

% Create table to fill values (adimensional unless stated in the title)
varNames=["file_name", "initial_velocity_cgs", "weber", "bond", "ohnesorge", "parent_folder", "number_of_harmonics", ...
    "max_width", "contact_time_ms", "coef_restitution", "north_pole_min_height", "north_pole_exp_min_height", ...
    "max_contact_radius", "spread_time_ms", "spread_time_width_ms", ...
    "contact_time_exp_ms", "coef_rest_exp", "max_contact_radius_exp", "spread_time_exp_ms", ...
    "energy_modes"];
varTypes=["string", "double", "double", "double", "double", "string", "double", "double", "double", "double", ...
    "double", "double", "double", "double", "double", "double", "double", "double", "double", 'cell'];

sz = [length(files_folder) length(varNames)];
fname = sprintf("postprocessing%s.mat", prefix);
if isfile(fullfile(root_folder, "2_output", fname))
    load(fullfile(root_folder,"2_output", fname), "data");
    
    
    data{length(files_folder)+1, 'file_name'} = missing;
    data = flipud(data);
    % Filter values which are not present in data
    files_folder = files_folder(cellfun(@(X) ~contains(X, data.file_name), {files_folder.name}));
else
    data = table();
end

if size(data, 2) ~= length(varNames) || sum(data.Properties.VariableNames ~= varNames) > 2
    data = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);
end

data = data(~isnan(data.initial_velocity_cgs), :);
%pixel = 5e-4; %Threshold for experimental contact
fnames = data{:, 1};
T = length(files_folder); 
if debug == true
    energies = table('Size', [1, 6], 'VariableTypes', ["double", "double", "double", "double", "double", "double"], ...
    'VariableNames', ["Weber", "KEin", "KEout", "PEout", "Eout", "epsilon"]);
end

% Sort by age (numeric field)
for ii = 1:length(files_folder)
    try
        % If it has been calculated before, or is it the postprocessing.mat
        % file, or it has eerror in its name, skip it. 
        if ismember(files_folder(ii).name, fnames) || ...
                contains(files_folder(ii).name, "postprocessing") || ...
                contains(lower(files_folder(ii).name), "error"); continue; 
        end
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

        if contains(lastwarn, 'not found'); error(lastwarn);  end
        
        % Extract dimensional variables sigma, rho, R, We
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
        
        pixel = 0.02 * length_unit; pixel_adim = pixel/length_unit;
        max_width = -inf;
        contact_time = nan; touch_time = nan; liftoff_time = nan;
        contact_time_exp = nan; liftoff_time_exp = nan;
        coef_restitution = nan; Vin = nan; Vout = nan;
        coef_restitution_exp = nan; Vout_exp = nan;
        north_pole_min_height = inf; north_pole_exp_min_height = inf;
        max_contact_radius = -inf; spread_time = nan; spread_time_width = nan;
        max_contact_radius_exp = -inf; spread_time_exp = nan; 
        
        velocity_unit = length_unit/time_unit;
        energy_constant = sigmaS/(rhoS * Ro * velocity_unit.^2);
        idxs = 1:PROBLEM_CONSTANTS.nb_harmonics;
        Xl = (2*pi./(idxs .* (2 * idxs + 1))); Yl = (2*pi * (idxs.^2 + idxs - 2)./(2*idxs+1));
        deformation_modes_energies = nan;
        A = (velocity_unit.^2/default_physical.initial_velocity.^2)/(2*pi/3);

        Vn = abs(default_physical.initial_velocity);
        g = default_physical.g;
        % Calculating experimental start of contact to substract (negative value)
        if g == 0
            t0 = 0.02*Ro/Vn;
        else
            t0 = real((Vn - sqrt(Vn^2 - 2*g*pixel_adim*Ro))/g); 
            if Vn^2 < 2*g*pixel_adim*Ro 
                fprintf("SKIPPING postprocessing (We = %.2e, Oh= %.2e, Bo = %.2e), Simul %d/%d \n", ...
                   Westar, Oh, Bo, ii, T);
                continue; 
            end            
        end
        V0 = sqrt(Vn^2-2*g*Ro*pixel_adim);
        X = @(t) pixel_adim*Ro - t*V0 -g*t^2/2;
        assert(abs(X(t0)) <= 1e-13 ); assert(abs(X(0) - pixel_adim*Ro) < 1e-13);
        fprintf("Starting postprocessing (We = %.2e, Oh= %.3e, Bo = %.3e), Simul %d/%d\n", ...
            Westar*(V0/Vn)^2, Oh, Bo, ii, T);
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
                    % Numerical variables
                    touch_time = recorded_times(jj);
                    Vin = recorded_conditions{jj}.center_of_mass_velocity;
                    CM_in = recorded_conditions{jj}.center_of_mass;
                    Ein = 1/2 * Vin^2;
                    % Experimental variables
                    Vin_exp = Vin + t0*g; assert(abs(Vin_exp+V0) < 1e-12);
                    %fprintf("Westar experimental is %.2e \n\n", Westar * Vin_exp^2/Vn^2);
                    Ein_exp = 1/2 * Vin_exp^2; % Now we are measuring 
                    touch_time_exp = touch_time - 1000*t0; % We shift time in miliseconds
                    CM_in_exp = (1+pixel_adim)* CM_in;
                    if debug == true 
                        energies(ii, 1:2) = {Westar * (V0/Vn)^2, Ein_exp};
                        fprintf("Input:  KE = %.3e.\n", Ein_exp); 
                    end
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
            north_pole_min_height = drop_height(0);
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
            
            north_pole_exp_min_height = min(north_pole_exp_min_height, north_pole_exp_max_projection);

            % Experimental contact_radius && coef restitution
            if isnan(liftoff_time_exp)
                % If contact ended numerically but simulation ended, record lift off time anyways
                if drop_height(pi) >= pixel/length_unit %||((size(recorded_conditions, 1)-1 == jj && ~isnan(liftoff_time)))
                    liftoff_time_exp = recorded_times(jj);
                    Vout_exp = recorded_conditions{jj}.center_of_mass_velocity;
                    Eout_exp = 1/2*Vout_exp^2 + (recorded_conditions{jj}.center_of_mass ...
                        - CM_in_exp)*g;

                    A = (velocity_unit.^2/default_physical.initial_velocity.^2)/(2*pi/3);
                    deformation_modes_energies = Xl .* (recorded_conditions{jj}.deformation_velocities(idxs)/velocity_unit).^2 + ...
                        energy_constant * Yl .* (recorded_conditions{jj}.deformation_amplitudes(idxs)/length_unit).^2;
                    deformation_modes_energies = A*deformation_modes_energies;
                    if debug == true
                        energies(ii, 3:6) = { 1/2*Vout_exp.^2, Eout_exp-1/2*Vout_exp.^2, Eout_exp, sqrt(abs(Eout_exp/Ein_exp))};
                        fprintf("Output: KE = %.3e, PE = %.3e, Eout = %.3e\n", 1/2*Vout_exp.^2, Eout_exp-1/2*Vout_exp.^2, Eout_exp); 
                    end
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
                contact_time_exp = liftoff_time_exp - touch_time_exp;
                coef_restitution_exp = sqrt(abs(Eout_exp/Ein_exp));
                
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
            kk = 1;
            theta2 = theta_vector(recorded_conditions{jj}.contact_points+kk);
            while drop_height(theta2) < pixel/length_unit
                kk = kk + 1;
                theta2 = theta_vector(recorded_conditions{jj}.contact_points + kk);
            end
            theta1 = pi;
            while abs(theta1-theta2) > pi/1e+4
                theta_middle = (theta1+theta2)/2;
                if drop_height(theta_middle) < pixel/length_unit
                    theta1 = theta_middle;
                else
                    theta2 = theta_middle;
                end
            end
            current_contact_radius_exp = sin(theta1) ...
                * drop_radius(theta1);
            
            if current_contact_radius_exp > max_contact_radius_exp
                max_contact_radius_exp = current_contact_radius_exp;
                spread_time_exp = recorded_times(jj) - touch_time_exp;
            end
        end

        %V0 = abs(default_physical.initial_velocity);
        %g = default_physical.g;
        %t0 = (-V0 + sqrt(V0^2 - 2*g*pixel_adim*Ro))/g; % Calculating experimental start of contact to substract
        %X = @(t) -t * V0 - g/2*t.^2 + (1 + pixel_adim)*Ro; 
        %rc = @(t) 10* sqrt(Ro^2 - (X(t) - 0.02*Ro).^2); % Contat radius in cm (10 because we go from cm to mm)
    
    
        % Add value to tables
        [~, dropLiquid, lol] = fileparts(files_folder(ii).folder); dropLiquid = strcat(dropLiquid, lol);
        data(ii, :) = {files_folder(ii).name, Vin_exp, Westar*Vin_exp^2/Vn^2, Bo, Oh, dropLiquid, PROBLEM_CONSTANTS.nb_harmonics, ...
            max_width, contact_time, coef_restitution, north_pole_min_height, ...
            north_pole_exp_min_height, max_contact_radius, spread_time, spread_time_width, ...
            contact_time_exp, coef_restitution_exp, max_contact_radius_exp, spread_time_exp, {deformation_modes_energies}};
        
    catch me
        if contains(files_folder(ii).name, "error"); continue; end
        fprintf("Exception at: %s %s \n", files_folder(ii).folder, files_folder(ii).name)
        disp(me);
        switch me.identifier
            case 'MATLAB:load:unableToReadMatFile'
                disp(" deletingfiles!");
                delete(fullfile(files_folder(ii).folder, files_folder(ii).name));
        end
    end
end
delete(gcp("nocreate")); % Deleting current parallel workers
% Filtering
data = rmmissing(data, 'DataVariables','file_name'); %data = data(data.max_width > -inf, :);

writetable(data, fullfile(root_folder, "2_output", replace(fname, ".mat", ".csv")));
s = fullfile(root_folder, "2_output", fname);
warning ('on','all');
if ~exist(s, "file")
    save(s, "data");
else
    warning('Did not save the postprocessing file. There is one in the directory already');
end
system('python3 sending_email.py');
diary off
% 
% hold off; scatter(energies.Weber, energies.KEout, 50, 'filled', 'DisplayName', 'KEout'); hold on
% scatter(energies.Weber, energies.PEout, 50, 'filled', 'DisplayName', 'PEout')
% scatter(energies.Weber, energies.Eout, 50, 'filled', 'DisplayName', 'KEout+PEout')
% yline(0); legend('show', 'FontSize', 14); xlabel("Weber");


function plotter(prefixString)

    % Get a list of all CSV files in the directory
    files = dir(fullfile("../../**/*.csv"));
    
    % Loop through each file
    for i = 1:length(files)
        % Get the file name
        fileName = files(i).name;
        
        % Check if the file name contains the prefix string
        if contains(fileName, prefixString)
            % Read the CSV file (assuming no headers)
            data = sortrows(readmatrix(fullfile(directory, fileName)));
            
            % Plot the data
            %figure;
            plot(data(:, 1), data(:, 2));  % Adjust columns as needed
            %title(['Plot of ', fileName]);
            %xlabel('X-axis');
            %ylabel('Y-axis');
        end
    end
end