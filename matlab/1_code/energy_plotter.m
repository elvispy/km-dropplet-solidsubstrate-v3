% Plots the deformation amplitudes, velocities and energies evolution for 
% a specific file to be chosen

function energy_plotter(varargin)
    % Add functions to calculate maximum width
    close all;
    %hold on;
    safe_folder = fullfile(fileparts(mfilename('fullpath')), "simulation_code");
    addpath(safe_folder, '-begin');
    filepathh = fullfile(fileparts(pwd), "2_output");
    [file, path] = uigetfile(filepathh, "*.m");
    fullfilepath = fullfile(path, file);
    %clear is_adim
    load(fullfilepath, "recorded_conditions", "default_physical");
    load(fullfilepath, "time_unit", "length_unit");
    postprocessingfile = fullfile(path, '..', 'postprocessing.mat');
    load(postprocessingfile, "data");

    %time_unit
    %velocity_unit = length_unit/time_unit;
    nb_harmonics = recorded_conditions{1}.nb_harmonics;
    %pressure_unit = default_physical.rhoS * default_physical.initial_velocity.^2;
    
    deformation_amplitudes = cell2mat(cellfun(@(x) x.deformation_amplitudes, ...
        recorded_conditions, 'UniformOutput',false))';
    deformation_velocities = cell2mat(cellfun(@(x) x.deformation_velocities, ...
        recorded_conditions, 'UniformOutput',false))';
    
    %Dimensional parameters
    COM = cell2mat(cellfun(@(x) x.center_of_mass, ...
        recorded_conditions, 'UniformOutput',false))';
    
    COM_velocities = cell2mat(cellfun(@(x) x.center_of_mass_velocity, ...
        recorded_conditions, 'UniformOutput',false))';

    times = cellfun(@(x) x.current_time, recorded_conditions)';

    
    cidx = find(cellfun(@(x) x.contact_points, recorded_conditions), 1, 'last'); 
    tidx = min(length(times), floor(1.5*cidx));
    times = times(1:tidx);
    deformation_velocities = deformation_velocities(:, 1:tidx);
    deformation_amplitudes = deformation_amplitudes(:, 1:tidx);
    
    sigmaS = default_physical.sigmaS; rhoS = default_physical.rhoS; Ro = default_physical.undisturbed_radius;
    g = default_physical.g;
    Westar = default_physical.rhoS * default_physical.initial_velocity^2 * ...
        default_physical.undisturbed_radius / default_physical.sigmaS;
    Oh = default_physical.nu * sqrt(default_physical.rhoS / (default_physical.sigmaS ...
        * default_physical.undisturbed_radius));
    Bo = rhoS * default_physical.g * Ro^2 / sigmaS;


    %% Plotting energy contributions
    energy_unit = 4*pi* Ro^2 * sigmaS;    
    
    %C = sigmaS * Ro^2/energy_unit;
    saving_figure_energy = figure('Position', [100, 500, 1900, 500]); % Wider than tall
    hold on;
    %recorded_conditions{ii}.amplitude_defor = 1;
    idxs = 1:nb_harmonics;

    Xl = (2*pi./(idxs .* (2 * idxs + 1)))'; Yl = (2*pi * (idxs.^2 + idxs - 2)./(2*idxs+1))';
    
    deformation_energies = rhoS * Ro^3 * Xl .* deformation_velocities(idxs, :).^2 + ...
        sigmaS * Yl .* deformation_amplitudes(idxs, :).^2;
    
    mechanical_energy = (4*pi/3 * rhoS * Ro^3 *(g*COM  + 1/2 * COM_velocities.^2));
    total_energy = (sum(deformation_energies(idxs, :), 1)+mechanical_energy);
    disp(any(total_energy(1:(end-1)) < total_energy(2:end)));

    plot(times/time_unit, sum(deformation_energies(idxs, :), 1)/energy_unit, 'LineWidth',2);
    plot(times/time_unit, mechanical_energy/energy_unit, 'LineWidth', 3);
    plot(times/time_unit, 4*pi/3 * rhoS * Ro^3 *(g*COM)/energy_unit, 'LineWidth', 3);
    plot(times/time_unit, (sum(deformation_energies(idxs, :), 1)+mechanical_energy)/energy_unit, 'LineWidth', 3);
    xline(times(cidx)/time_unit, 'LineWidth', 2);
    grid on;
    % Only show a subset of 
    set(gca, 'FontSize', 16);
    legend(["Total Surface energy", "Total Mechanical energy COM", "Potential energy only", "Surface + Mechanical", "Contact ends"], 'FontSize', 18, "Location", "SouthEastOutside");
    xlabel('$ t/t_\sigma $', 'Interpreter','Latex', 'FontSize', 20);
    ylabel('$ E/(4 \pi R_o^2 \sigma) $', 'Interpreter','Latex', 'FontSize', 20);
   
    title(sprintf("Energy contribution per mode with We = %.3g, Oh = %.3g, Bo = %.3g", Westar, Oh, Bo));
    yl = get(gca, 'YLim'); yl = [-.01, max(abs(yl))];
    %set(gca, 'YLim', yl); %set(gca, 'XLim', 1.05*get(gca, 'XLim'))
    saveas(saving_figure_energy, fullfile(path, replace(file, ".mat", "_energy.png")));

    
end




