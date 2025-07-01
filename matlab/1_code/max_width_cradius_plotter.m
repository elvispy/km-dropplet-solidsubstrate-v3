% This script will plot a normal animation, but with the evolution
% of the contact radius, max width, and minimum experimental height
function max_width_cradius_plotter(varargin)
    if nargin == 1
       
    else
        % Add functions to calculate maximum width
        safe_folder = fullfile(fileparts(mfilename('fullpath')), "simulation_code");
        addpath(safe_folder, '-begin');
        [file, path] = uigetfile("../2_output/*.mat", 'MultiSelect', 'on');
        if ischar(file) || isstr(file); file = {file}; end
        fullfilepath = cellfun(@(file) fullfile(path, file), file, 'UniformOutput',false);
        physicalParameters = cell(length(file), 1);
        units = cell(length(file), 1);
        
        for ii = 1:length(file)
            physicalParameters{ii} = load(fullfilepath{ii}, ...
                "recorded_conditions", "default_physical", "PROBLEM_CONSTANTS", "recorded_times");
            units{ii} = load(fullfilepath{ii}, "length_unit", "velocity_unit", "pressure_unit"); 
            vidObj = replace(fullfilepath{ii}, ".mat", ".mp4");
        end
        vidObj = VideoWriter(vidObj, "MPEG-4");
        set(vidObj, 'Quality', 100, 'FrameRate', 10);
        open(vidObj);
        close all;
        pplots = figure(1); set(pplots, 'Position', [954 94 900 600]);
        
        
        texts = cell(length(file), 1);
        max_widths_tracker = ones(length(file), 1); spread_time_tracker = zeros(length(file), 1);
        north_pole_exp_min = 3* ones(length(file), 1);
        max_contact_radius_tracker = zeros(length(file), 1); spread_widths_tracker = zeros(length(file), 1);
        theta_vector = physicalParameters{1}.PROBLEM_CONSTANTS.theta_vector;
        time_unit = arrayfun(@(jj) sqrt(physicalParameters{jj}.default_physical.rhoS* physicalParameters{jj}.default_physical.undisturbed_radius^3/physicalParameters{jj}.default_physical.sigmaS), ...
            1:length(file), 'UniformOutput', false);
        times = arrayfun(@(jj) physicalParameters{jj}.recorded_times, ...
            1:length(file), 'UniformOutput', false);
        sz = min(arrayfun(@(jj) size(physicalParameters{jj}.recorded_times, 1), 1:length(file)));
        n_shots = min(200, sz); nbefore = 30;
        max_widths  = 2*ones(length(file), n_shots); contact_radii02R = zeros(length(file), n_shots); cradiibefore = zeros(length(file), nbefore);
        north_poles = 2*ones(length(file), n_shots); contact_radii00R = zeros(length(file), n_shots);
        max_radii_time = zeros(length(file), 1); max_contact_radii = zeros(length(file), 1);
        times_idx = ceil(linspace(1, sz, n_shots));
        wes = zeros(length(file), 1); ohs = zeros(length(file), 1); bos = zeros(length(file), 1);
        for idx = 1:n_shots
            ii = times_idx(idx);
            for jj = 1:length(physicalParameters)
                %hold on;
                adim_conditions = physicalParameters{jj}.recorded_conditions{ii};
                %default_physical = physicalParameters{jj}.default_physical;
                
                adim_conditions.center_of_mass = adim_conditions.center_of_mass/units{jj}.length_unit;
                adim_conditions.pressure_amplitudes = adim_conditions.pressure_amplitudes/units{jj}.pressure_unit;
                adim_conditions.center_of_mass_velocity = adim_conditions.center_of_mass_velocity/units{jj}.velocity_unit;
                adim_conditions.deformation_amplitudes = adim_conditions.deformation_amplitudes/units{jj}.length_unit;
                
                
                adim_deformations = adim_conditions.deformation_amplitudes;
                adim_CM = adim_conditions.center_of_mass;
                drop_radius = zeta_generator(adim_deformations);
                drop_radius = @(theta) 1 + drop_radius(theta);
                drop_height = @(theta) cos(theta) .* drop_radius(theta) + adim_CM;
                %tic = sqrt(default_physical.rhoS* default_physical.undisturbed_radius^3/default_physical.sigmaS);
                
                %recorded_conditions{ii}.amplitude_defor = 1;
                pixel_adim = 0.02;
                default_physical = physicalParameters{jj}.default_physical;
                if idx <= 2
                    
                    g = default_physical.g; Ro = default_physical.undisturbed_radius; 
                    Vn = abs(default_physical.initial_velocity);
                    % Calculating experimental start of contact to substract (negative value)
                    if g == 0
                        t0 = 0.02*Ro/Vn;
                    else
                        t0 = real((Vn - sqrt(Vn^2 - 2*g*pixel_adim*Ro))/g); 
                        if Vn^2 < 2*g*pixel_adim*Ro 
                            fprintf("Invalid data");
                            error("Not a valid impact velocity. Increase it.");
                        end            
                    end
                    V0 = sqrt(Vn^2-2*g*Ro*pixel_adim);
                    X = @(t) pixel_adim*Ro - t*V0 - g*t.^2/2;
                    tt0 = linspace(-t0, 0, nbefore);
                    cradiibefore(jj, :) = sin(acos((X(tt0 + t0) + (1-pixel_adim)*Ro)./Ro)); % Non dimensional contact radius!
                
                
                
                    wes(jj) = default_physical.rhoS * V0^2 * ...
                        default_physical.undisturbed_radius / default_physical.sigmaS;
                    ohs(jj) = default_physical.nu * sqrt(default_physical.rhoS / (default_physical.sigmaS ...
                        * default_physical.undisturbed_radius));
                    bos(jj) = default_physical.rhoS * default_physical.g * default_physical.undisturbed_radius^2 / default_physical.sigmaS;

                end
                
                %% Calculating min height (experimental side-view projection)
                north_pole_exp_min_height = drop_height(0); currang = pi/4; dtheta = pi/2; N = 9;
                while dtheta > pi/N^3
                   values = drop_height(linspace(currang-dtheta/2, currang+dtheta/2, N));
                   [maxval, maxindex] = max(values);
                   if maxval > north_pole_exp_min_height
                       north_pole_exp_min_height = maxval;
                   end
                   currang = max(currang - dtheta/2 + (maxindex-1)/N * dtheta, dtheta/(2*N));
                   dtheta = dtheta/N;
                   if north_pole_exp_min_height == drop_height(0)
                       break;
                   end
                end
                
                north_poles(jj, idx) = north_pole_exp_min_height;
                north_pole_exp_min(jj) = min(north_pole_exp_min(jj), north_pole_exp_min_height);
                

                %% Calculating (current) contact radius
                % 0.02 Treshold
                %pixel_adim = 0.02;
                
                kk = 1;
                current_contact_points_exp = adim_conditions.contact_points;
                theta2 = theta_vector(current_contact_points_exp +kk);
                while drop_height(theta2) < pixel_adim
                    kk = kk + 1;
                    theta2 = theta_vector(current_contact_points_exp + kk);
                end
                theta1 = pi;
                while abs(theta1-theta2) > pi/1e+4
                    theta_middle = (theta1+theta2)/2;
                    if drop_height(theta_middle) < pixel_adim
                        theta1 = theta_middle;
                    else
                        theta2 = theta_middle;
                    end
                end
                current_contact_radius_exp = sin(theta1) ...
                    * drop_radius(theta1);
                contact_radii02R(jj, idx) = current_contact_radius_exp;
                if current_contact_radius_exp > max_contact_radius_tracker(jj)
                    max_contact_radius_tracker(jj) = current_contact_radius_exp;
                    spread_time_tracker(jj) = times{jj}(ii);
                end
                
                % 0.00R treshhold
                pixel_adim = 1e-4;
                
                kk = 1;
                current_contact_points_exp = adim_conditions.contact_points;
                theta2 = theta_vector(current_contact_points_exp +kk);
                while drop_height(theta2) < pixel_adim
                    kk = kk + 1;
                    theta2 = theta_vector(current_contact_points_exp + kk);
                end
                theta1 = pi;
                while abs(theta1-theta2) > pi/1e+4
                    theta_middle = (theta1+theta2)/2;
                    if drop_height(theta_middle) < pixel_adim
                        theta1 = theta_middle;
                    else
                        theta2 = theta_middle;
                    end
                end
                current_contact_radius_exp = sin(theta1) ...
                    * drop_radius(theta1);
                contact_radii00R(jj, idx) = current_contact_radius_exp;
                
%                 if current_contact_radius_exp > max_contact_radii(jj)
%                     max_radii_time(jj) = times{jj}(ii);
%                     max_contact_radii(jj) = current_contact_radius_exp;
%                 end
                
                %% Calculating (current) max width
                current_max_width_adim = maximum_contact_radius(adim_conditions.deformation_amplitudes);
                max_widths(jj, idx) = current_max_width_adim;
                if current_max_width_adim > max_widths_tracker(jj)
                    max_widths_tracker(jj) = current_max_width_adim;
                    spread_widths_tracker(jj) = times{jj}(ii);
                end
                
                if false
                    % Plotting
                    % Current time, maximum radius, maximum equatorial radius
                    ax = subplot(1, length(file), jj);
                    plot_condition(ax, adim_conditions, 1.75); %changed from jj==1 to true (now that are subplots)
                    texts{jj} = sprintf('$t = %.2f$  (ms),$ t_{r_c} = %.2f$  (ms), $ t_m = %.2f$ (ms)', 1000*times{jj}(ii), 1000* spread_time_tracker(jj), ...
                        1000* spread_widths_tracker(jj));


                    %% Extreme values plotter
                    yline(north_pole_exp_min(jj), 'r--', 'LineWidth', 4); % North pole min height
                    xline(-max_widths_tracker(jj), '--', 'LineWidth', 4); xline(max_widths_tracker(jj), '--', 'LineWidth', 4); % Max width
                    plot([-max_contact_radius_tracker(jj), max_contact_radius_tracker(jj)], [0, 0], ...
                        '-', 'LineWidth', 1+1.5*jj);

                    title(sprintf('We = %.2e',wes(jj)), 'FontSize', 14);

                    %% Plotting everything
                    set(gca, 'FontSize', 20);
                    fill([-10, 10, 10, -10], [0, 0, -1, -1] * 1e6, [235 176 0]./256, 'EdgeColor', 'none', 'FaceAlpha',0.3);
                    text(0.02, 0.95, texts{jj}, 'Interpreter', 'latex', 'Units', 'normalized', ...
                        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
                        'FontSize', 20, 'FontWeight', 'bold');
                    %title(sprintf('t = %.4f, CM velocity = %.2f (cm/s), Contact Points = %g',...
                    %    adim_conditions.current_time/tic, physicalParameters{ii}.center_of_mass_velocity, ...
                    %    adim_conditions.contact_points), ...
                    %    'FontSize', 24);
                    grid on;
                    set(gcf,'Color','white');
                    xlabel('$ x/R_o $', 'Interpreter','Latex', 'FontSize', 28);
                    ylabel('$ y/R_o $', 'Interpreter','Latex', 'FontSize', 28);
                end
            end % end inner for (videos)
            
            if false; writeVideo(vidObj, getframe(gcf)); end
            
        end % end outer for (video
        close(vidObj);
        figure(2); set(gcf, 'Position', [477 595 1572 420]);
        for jj = 1:length(physicalParameters)
            subplot(1, 3, 1); hold on; grid on; set(gca, 'FontSize', 16);
            title("Minimum Height", 'FontSize', 16); xlabel("Time (ms)", 'FontSize', 14);
            plot(1000*times{jj}(times_idx), north_poles(jj, :), '--'    , 'LineWidth', jj+1, 'DisplayName', sprintf("We = %.2e, Oh = %.1e", wes(jj), ohs(jj)));
            legend('FontSize', 14);
            subplot(1, 3, 2); hold on; grid on; set(gca, 'FontSize', 16);
            title("Contact Radius", 'FontSize', 16); xlabel("Time (ms)", 'FontSize', 14);
            plot(1000*tt0, cradiibefore(jj, :), 'LineWidth', jj+1, 'DisplayName', "");
            plot(1000*times{jj}(times_idx), contact_radii02R(jj, :), 'LineWidth', jj+1, 'DisplayName', sprintf("0.02R, We = %.2e, Oh = %.1e", wes(jj), ohs(jj)));
            plot(1000*times{jj}(times_idx), contact_radii00R(jj, :), 'LineWidth', jj+1, 'DisplayName', sprintf("0.00R, We = %.2e, Oh = %.1e", wes(jj), ohs(jj)));
            scatter(1000*spread_time_tracker(jj), max_contact_radius_tracker(jj), "*", 'DisplayName', sprintf('Time to maximum radius We = %.2e, Oh = %.1e', wes(jj), ohs(jj)));
            legend('FontSize', 12, 'Location', 'southeast');
            subplot(1, 3, 3); hold on; grid on; set(gca, 'FontSize', 16);
            title("Maximum Width", 'FontSize', 16); xlabel("Time (ms)", 'FontSize', 14);
            plot(1000*times{jj}(times_idx), max_widths(jj, :), '--', 'LineWidth', jj+1, 'DisplayName', sprintf("We = %.2e, Oh = %.1e", wes(jj), ohs(jj)));
            legend('FontSize', 14);
            
        end
        
        % Create a new Excel file

        % Create a new Excel file
        excelFileName = '../2_output/Simulation_Results.xlsx';

        % Loop through each simulation file
        for jj = 1:length(file)
            % Use the full file name (including .mat) as the sheet name
            sheetName = file{jj}; % This will include the .mat extension

            % Prepare the data for the current sheet
            
            timeData = 1000 .* [t0+tt0'; t0+times{jj}(times_idx)]; % Convert time to milliseconds
            Ro = physicalParameters{jj}.default_physical.undisturbed_radius; % Get Ro for dimensional conversion
            contactRadiiData = Ro/100 .* [cradiibefore(jj, :)'; contact_radii02R(jj, :)']; % Convert to meters (dimensional)

            % Combine time and contact radii data into a table
            dataTable = table(timeData, contactRadiiData, 'VariableNames', {'Time_ms', 'Contact_Radii_m'});

            % Write the time and contact radii data to the Excel sheet
            writetable(dataTable, excelFileName, 'Sheet', sheetName, 'Range', 'A1');

            % Create a small table for Weber, Ohnesorge, and Bond numbers
            weOhBoTable = table({'Weber'; 'Ohnesorge'; 'Bond'}, [wes(jj); ohs(jj); bos(jj)], ...
                'VariableNames', {'Parameter', 'Value'});

            % Write the Weber, Ohnesorge, and Bond table to the Excel sheet
            writetable(weOhBoTable, excelFileName, 'Sheet', sheetName, 'Range', 'D1');
        end

        disp(['Data exported to ' excelFileName]);

        
    end % end nargin == 1

end % end function definition