% Plots the deformation amplitudes, velocities and energies evolution for 
% a specific file to be chosen

function energy_plotter(varargin)
    % Add functions to calculate maximum width
    close all;
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
    %MM = 200; 
    %H = floor(linspace(1, recorded_conditions{1}.nb_harmonics, 20));
    %pressure_amplitudes = zeros(length(H), MM);
    %deformation_amplitudes = zeros(length(H), MM);
    %times = zeros(MM, 1);
    velocity_unit = length_unit/time_unit;
    nb_harmonics = recorded_conditions{1}.nb_harmonics;
    pressure_unit = default_physical.rhoS * default_physical.initial_velocity.^2;
    pressure_amplitudes = cell2mat(cellfun(@(x) x.pressure_amplitudes, ...
        recorded_conditions, 'UniformOutput',false))'/pressure_unit;
    deformation_amplitudes = cell2mat(cellfun(@(x) x.deformation_amplitudes, ...
        recorded_conditions, 'UniformOutput',false))'/length_unit;
    deformation_velocities = cell2mat(cellfun(@(x) x.deformation_velocities, ...
        recorded_conditions, 'UniformOutput',false))'/(length_unit/time_unit);
    
    %Dimensional parameters
    COM = cell2mat(cellfun(@(x) x.center_of_mass, ...
        recorded_conditions, 'UniformOutput',false))';
    
    COM_velocities = cell2mat(cellfun(@(x) x.center_of_mass_velocity, ...
        recorded_conditions, 'UniformOutput',false))';

    times = cellfun(@(x) x.current_time, recorded_conditions)'/time_unit;
    % To find maximum contact points
    %[~, cidx] = max(cellfun(@(x) x.contact_points, recorded_conditions)); 
    % To find contact time
    cidx = find(cellfun(@(x) x.contact_points, recorded_conditions), 1, 'last'); 
    tidx = min(length(times), floor(1.5*cidx));
    times = times(1:tidx);
    deformation_velocities = deformation_velocities(:, 1:tidx);
    deformation_amplitudes = deformation_amplitudes(:, 1:tidx);
    pressure_amplitudes    = pressure_amplitudes(:, 1:tidx);
    
    sigmaS = default_physical.sigmaS; rhoS = default_physical.rhoS; Ro = default_physical.undisturbed_radius;
    Westar = default_physical.rhoS * default_physical.initial_velocity^2 * ...
        default_physical.undisturbed_radius / default_physical.sigmaS;
    Oh = default_physical.nu * sqrt(default_physical.rhoS / (default_physical.sigmaS ...
        * default_physical.undisturbed_radius));
    Bo = rhoS * default_physical.g * Ro^2 / sigmaS;
    
    % Plotting amplitudes
    if false
    saving_figure = figure('Position', [100, 300, 1000, 300]); % Wider than tall
    hold on;
    %recorded_conditions{ii}.amplitude_defor = 1;
    idxs = 1:15;
    idxs = idxs(max(abs(deformation_amplitudes(1:idxs(end), :)), [] , 2) > 5e-3); %2.^(1:floor(log2(nb_harmonics)));
    %cmap = colormap('spring'); %disp(size(cmap));
    %numColors = size(cmap, 1);
    lol = jet(length(idxs)+2);
    colororder(flipud(lol(3:end, :)));%cmap(floor(linspace(1, numColors, length(idxs))), :));
    %disp(size(times)); disp(size(deformation_amplitudes));
    plot(times, (deformation_amplitudes(idxs, :)), 'LineWidth',2);
    %plot(times, (deformation_velocities(idxs, :)), 'LineWidth',2);
    %xline(times(cidx), 'LineWidth', 2);
    % Only show a subset of 
    set(gca, 'FontSize', 16);
    legend(arrayfun(@(i) string(i), idxs), 'FontSize', 12);
    %grid on;
    xlabel('$ t/t_s $', 'Interpreter','Latex', 'FontSize', 20);
    ylabel('$ r/R_o $', 'Interpreter','Latex', 'FontSize', 20);
    
    
    
    title(sprintf("Deformation Amplitudes with We = %.3g, Oh = %.3g, Bo = %.3g", Westar, Oh, Bo));
    yl = get(gca, 'YLim'); yl = [-max(abs(yl)), max(abs(yl))];
    set(gca, 'YLim', yl); set(gca, 'XLim', 1.05*get(gca, 'XLim'))
    end
    
    

    %% Plotting pressures
    if false
    saving_figure_pressure = figure('Position', [100, 50, 1000, 300]); % Wider than tall
    hold on;
    %recorded_conditions{ii}.amplitude_defor = 1;
    idxs = 1:(nb_harmonics+1);
    %idxs = idxs(max(abs(pressure_amplitudes), [] , 2) < 5); %2.^(1:floor(log2(nb_harmonics)));
    %cmap = colormap('spring'); %disp(size(cmap));
    %numColors = size(cmap, 1);
    lol = jet(length(idxs)+2);
    colororder(flipud(lol(2:end, :)));%cmap(floor(linspace(1, numColors, length(idxs))), :));
    %disp(size(times)); disp(size(deformation_amplitudes));
    plot(times, pressure_amplitudes(idxs(1:(end-1)), :), 'LineWidth',2);
    % Only show a subset of 
    set(gca, 'FontSize', 16);
    legend(arrayfun(@(i) string(i), idxs), 'FontSize', 12);
    xlabel('$ t/t_s $', 'Interpreter','Latex', 'FontSize', 20);
    ylabel('$ \frac{p}{\rho V_o^2} $', 'Interpreter','Latex', 'FontSize', 20, 'Rotation',90);
    
    sigmaS = default_physical.sigmaS; rhoS = default_physical.rhoS; Ro = default_physical.undisturbed_radius;
    Westar = default_physical.rhoS * default_physical.initial_velocity^2 * ...
        default_physical.undisturbed_radius / default_physical.sigmaS;
    Oh = default_physical.nu * sqrt(default_physical.rhoS / (default_physical.sigmaS ...
        * default_physical.undisturbed_radius));
    Bo = rhoS * default_physical.g * Ro^2 / sigmaS;

    title(sprintf("Deformation Amplitudes with We = %.3g, Oh = %.3g, Bo = %.3g", Westar, Oh, Bo));
    yl = get(gca, 'YLim'); yl = [-max(abs(yl)), max(abs(yl))];
    set(gca, 'YLim', yl); set(gca, 'XLim', 1.05*get(gca, 'XLim'))
    
    end
    %cd(curr);
    %saveas(saving_figure, "../../0_data/manual/pressure_plotter", 'fig');
    %print(saving_figure, '-depsc', '-r300', "../../0_data/manual/pressure_plotter.eps");


    %% Plotting energy contributions
    energy_unit = rhoS * Ro^3 * velocity_unit^2;
    
    C = sigmaS * Ro^2/energy_unit;
    saving_figure_energy = figure('Position', [100, 500, 1000, 300]); % Wider than tall
    hold on;
    %recorded_conditions{ii}.amplitude_defor = 1;
    idxs = 1:nb_harmonics;
    %idxs = idxs(max(abs(deformation_amplitudes), [] , 2) > 5e-3); %2.^(1:floor(log2(nb_harmonics)));
    %cmap = colormap('spring'); %disp(size(cmap));
    %numColors = size(cmap, 1);
    %lol = jet(length(idxs)+2);
    %colororder(flipud(lol(3:end, :)));%cmap(floor(linspace(1, numColors, length(idxs))), :));
    %disp(size(times)); disp(size(deformation_amplitudes));
    Xl = (2*pi./(idxs .* (2 * idxs + 1)))'; Yl = (2*pi * (idxs.^2 + idxs - 2)./(2*idxs+1))';
    deformation_energies = Xl .* deformation_velocities(idxs, :).^2 + ...
        C * Yl .* deformation_amplitudes(idxs, :).^2;
    
    mechanical_energy = (rhoS * Ro^3 *(COM - 0.02*Ro + 1/2 * COM_velocities.^2))/(energy_unit);
    A = (velocity_unit.^2/default_physical.initial_velocity.^2)/(2*pi/3);
    plot(times, A*sum(deformation_energies(idxs, :), 1), 'LineWidth',2);
    plot(times, mechanical_energy, 'LineWidth', 3);
    xline(times(cidx), 'LineWidth', 2);
    grid on;
    % Only show a subset of 
    set(gca, 'FontSize', 16);
    legend(["Total Surface energy", "Total Mechanical energy COM", "Contact ends"], 'FontSize', 12);
    xlabel('$ t/t_s $', 'Interpreter','Latex', 'FontSize', 20);
    ylabel('$ E/(\rho R_o^3 V^2) $', 'Interpreter','Latex', 'FontSize', 20);
   
    title(sprintf("Energy contribution per mode with We = %.3g, Oh = %.3g, Bo = %.3g", Westar, Oh, Bo));
    yl = get(gca, 'YLim'); yl = [-.1, max(abs(yl))];
    set(gca, 'YLim', yl); %set(gca, 'XLim', 1.05*get(gca, 'XLim'))
    saveas(saving_figure_energy, fullfile(path, replace(file, ".mat", "_energy.png")));

    
end

% function load_vars(str)
%     global errored
% 
%     if errored == true
%         str = "errored_" + str; 
%     end
% 
%     vars = load(str);
%     fn = fieldnames(vars);
%     for ii = 1:length(fn)
%         assignin('caller', fn{ii}, vars.(fn{ii}));
%     end
% 
% end



