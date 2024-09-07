function export_animation(varargin)
    if nargin == 1
       
    else
        % Add functions to calculate maximum width
        safe_folder = fullfile(fileparts(mfilename('fullpath')), "simulation_code");
        addpath(safe_folder, '-begin');
        [file, path] = uigetfile("*.mat");
        fullfilepath = fullfile(path, file);
        clear is_adim
        load(fullfilepath, "recorded_conditions", "is_adim");
        
        vidObj = VideoWriter(replace(fullfilepath, ".mat", ".mp4"), "MPEG-4");
        set(vidObj, 'Quality', 100, 'FrameRate', 30);
        open(vidObj);
        
        for ii = 1:size(recorded_conditions, 1)
            adim_conditions = recorded_conditions{ii};
            if exist('is_adim') == 0 || is_adim == false
                load(fullfilepath, "length_unit", "velocity_unit", "pressure_unit"); 
                adim_conditions.center_of_mass = adim_conditions.center_of_mass/length_unit;
                adim_conditions.pressure_amplitudes = adim_conditions.pressure_amplitudes/pressure_unit;
                adim_conditions.center_of_mass_velocity = adim_conditions.center_of_mass_velocity/velocity_unit;
                adim_conditions.deformation_amplitudes = adim_conditions.deformation_amplitudes/length_unit;
            end
            
            %recorded_conditions{ii}.amplitude_defor = 1;
            plot_condition(1, adim_conditions, 2);
            title(sprintf('t = %.2g, CP = %g',...
                adim_conditions.current_time, adim_conditions.contact_points))
            writeVideo(vidObj, getframe(gcf));
            
        end
        close(vidObj);
    end

end