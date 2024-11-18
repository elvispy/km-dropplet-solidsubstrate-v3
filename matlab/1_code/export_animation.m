function export_animation(varargin)
    if nargin == 1
       
    else
        % Add functions to calculate maximum width
        safe_folder = fullfile(fileparts(mfilename('fullpath')), "simulation_code");
        addpath(safe_folder, '-begin');
        [file, path] = uigetfile("*.mat");
        fullfilepath = fullfile(path, file);
        clear is_adim
        load(fullfilepath, "recorded_conditions", "is_adim", "default_physical");
        
        vidObj = VideoWriter(replace(fullfilepath, ".mat", ".mp4"), "MPEG-4");
        set(vidObj, 'Quality', 100, 'FrameRate', 40);
        open(vidObj);
        
        for ii = floor(linspace(1, size(recorded_conditions, 1), 600))
            adim_conditions = recorded_conditions{ii};
            if exist('is_adim') == 0 || is_adim == false
                load(fullfilepath, "length_unit", "velocity_unit", "pressure_unit"); 
                adim_conditions.center_of_mass = adim_conditions.center_of_mass/length_unit;
                adim_conditions.pressure_amplitudes = adim_conditions.pressure_amplitudes/pressure_unit;
                adim_conditions.center_of_mass_velocity = adim_conditions.center_of_mass_velocity/velocity_unit;
                adim_conditions.deformation_amplitudes = adim_conditions.deformation_amplitudes/length_unit;
            end
            tic = sqrt(default_physical.rhoS* default_physical.undisturbed_radius^3/default_physical.sigmaS);
            %recorded_conditions{ii}.amplitude_defor = 1;
            plot_condition(1, adim_conditions, 1.25);
            fill([-10, 10, 10, -10], [0, 0, -1, -1] * 1e6, [235 176 0]./256, 'EdgeColor', 'none', 'FaceAlpha',0.3);

            % title(sprintf('t = %.4f (ms), CM velocity = %.2f (cm/s), Contact Points = %g',...
            %     adim_conditions.current_time*1000, recorded_conditions{ii}.center_of_mass_velocity, ...
            %     adim_conditions.contact_points), ...
            %     'FontSize', 14);
            title(sprintf('t = %.4f, CM velocity = %.2f (cm/s), Contact Points = %g',...
                adim_conditions.current_time/tic, recorded_conditions{ii}.center_of_mass_velocity, ...
                adim_conditions.contact_points), ...
                'FontSize', 14);
            set(gca, 'FontSize', 16);
            xlabel('$ x/R_o $', 'Interpreter','Latex', 'FontSize', 20);
            ylabel('$ y/R_o $', 'Interpreter','Latex', 'FontSize', 20);
            writeVideo(vidObj, getframe(gcf));
            
        end
        close(vidObj);
    end

end