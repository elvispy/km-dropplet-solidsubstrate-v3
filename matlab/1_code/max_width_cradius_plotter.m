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
        
        texts = cell(length(file), 1);
        max_widths_tracker = ones(length(file), 1); spread_time_tracker = zeros(length(file), 1);
        max_contact_radius_tracker = zeros(length(file), 1); spread_widths_tracker = zeros(length(file), 1);
        theta_vector = physicalParameters{1}.PROBLEM_CONSTANTS.theta_vector;
        adim_times = arrayfun(@(jj) physicalParameters{jj}.recorded_times/sqrt(physicalParameters{jj}.default_physical.rhoS* physicalParameters{jj}.default_physical.undisturbed_radius^3/physicalParameters{jj}.default_physical.sigmaS), ...
            1:length(file), 'UniformOutput', false);
        sz = min(arrayfun(@(jj) size(physicalParameters{jj}.recorded_times, 1), 1:length(file)));
        for ii = ceil(linspace(1, sz, 200))
            
            for jj = 1:length(physicalParameters)
                %hold on;
                adim_conditions = physicalParameters{jj}.recorded_conditions{ii};
                %default_physical = physicalParameters{jj}.default_physical;
                
                adim_conditions.center_of_mass = adim_conditions.center_of_mass/units{jj}.length_unit;
                adim_conditions.pressure_amplitudes = adim_conditions.pressure_amplitudes/units{jj}.pressure_unit;
                adim_conditions.center_of_mass_velocity = adim_conditions.center_of_mass_velocity/units{jj}.velocity_unit;
                adim_conditions.deformation_amplitudes = adim_conditions.deformation_amplitudes/units{jj}.length_unit;
                pixel_adim = 0.02;
                %tic = sqrt(default_physical.rhoS* default_physical.undisturbed_radius^3/default_physical.sigmaS);
                
                %recorded_conditions{ii}.amplitude_defor = 1;
                plot_condition(1, adim_conditions, 1.45, jj==1);

                % Calculating (current) contact radius
                adim_CM = adim_conditions.center_of_mass;
                drop_radius = zeta_generator(adim_conditions.deformation_amplitudes);
                drop_radius = @(theta) 1 + drop_radius(theta);
                drop_height = @(theta) cos(theta) .* drop_radius(theta) + adim_CM;

                kk = 1;
                current_contact_points_exp = adim_conditions.contact_points;
                theta2 = theta_vector(current_contact_points_exp +kk);
                while drop_height(theta2) < pixel_adim
                    kk = kk + 1;
                    theta2 = theta_vector(current_contact_points_exp + kk);
                end
                theta1 = theta_vector(current_contact_points_exp + kk);
                while abs(theta1-theta2) > pi/1e+3
                    theta_middle = (theta1+theta2)/2;
                    if drop_height(theta_middle) < pixel_adim
                        theta1 = theta_middle;
                    else
                        theta2 = theta_middle;
                    end
                end
                current_contact_radius_exp = sin(theta1) ...
                    * drop_radius(theta1);
                if current_contact_radius_exp > max_contact_radius_tracker(jj)
                    max_contact_radius_tracker(jj) = current_contact_radius_exp;
                    spread_time_tracker(jj) = adim_times{jj}(ii);
                end

                %Calculating (current) max width
                current_max_width_adim = maximum_contact_radius(adim_conditions.deformation_amplitudes);
                if current_max_width_adim > max_widths_tracker(jj)
                    max_widths_tracker(jj) = current_max_width_adim;
                    spread_widths_tracker(jj) = adim_times{jj}(ii);
                end
                texts{jj} = sprintf('$t = %.2f, t_c = %.2f, t_m = %.2f$', adim_times{jj}(ii), spread_time_tracker(jj), ...
                    spread_widths_tracker(jj));

                xline(-max_widths_tracker(jj), '--', 'LineWidth', 4); xline(max_widths_tracker(jj), '--', 'LineWidth', 4);
                plot([-max_contact_radius_tracker(jj), max_contact_radius_tracker(jj)], [0, 0], ...
                    '-', 'LineWidth', 1+1.5*jj);

            end % end inner for (videos)
            
                % title(sprintf('t = %.4f (ms), CM velocity = %.2f (cm/s), Contact Points = %g',...
                %     adim_conditions.current_time*1000, recorded_conditions{ii}.center_of_mass_velocity, ...
                %     adim_conditions.contact_points), ...
                %     'FontSize', 14);
            
            set(gca, 'FontSize', 20);
            fill([-10, 10, 10, -10], [0, 0, -1, -1] * 1e6, [235 176 0]./256, 'EdgeColor', 'none', 'FaceAlpha',0.3);
            text(0.02, 0.95, texts, 'Interpreter', 'latex', 'Units', 'normalized', ...
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
            writeVideo(vidObj, getframe(gcf));
            
        end % end outer for (video
        close(vidObj);
    end %

end % end function definition