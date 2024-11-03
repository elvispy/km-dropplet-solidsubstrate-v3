function amplitudes_plotter(varargin)
    % Add functions to calculate maximum width
    close all;
    safe_folder = fullfile(fileparts(mfilename('fullpath')), "simulation_code");
    addpath(safe_folder, '-begin');
    [file, path] = uigetfile("*.mat");
    fullfilepath = fullfile(path, file);
    %clear is_adim
    load(fullfilepath, "recorded_conditions", "default_physical");
    load(fullfilepath, "time_unit", "length_unit");
    %MM = 200; 
    %H = floor(linspace(1, recorded_conditions{1}.nb_harmonics, 20));
    %pressure_amplitudes = zeros(length(H), MM);
    %deformation_amplitudes = zeros(length(H), MM);
    %times = zeros(MM, 1);
    nb_harmonics = recorded_conditions{1}.nb_harmonics;
    pressure_unit = default_physical.rhoS * default_physical.initial_velocity.^2;
    pressure_amplitudes = cell2mat(cellfun(@(x) x.pressure_amplitudes, ...
        recorded_conditions, 'UniformOutput',false))'/pressure_unit;
    deformation_amplitudes = cell2mat(cellfun(@(x) x.deformation_amplitudes, ...
        recorded_conditions, 'UniformOutput',false))'/length_unit;
    times = cellfun(@(x) x.current_time, recorded_conditions)'/time_unit;


    saving_figure = figure('Position', [100, 300, 1000, 300]); % Wider than tall
    hold on;
    %recorded_conditions{ii}.amplitude_defor = 1;
    idxs = 1:nb_harmonics;
    idxs = idxs(max(abs(deformation_amplitudes), [] , 2) > 5e-3); %2.^(1:floor(log2(nb_harmonics)));
    %cmap = colormap('spring'); %disp(size(cmap));
    %numColors = size(cmap, 1);
    lol = jet(length(idxs)+2);
    colororder(flipud(lol(3:end, :)));%cmap(floor(linspace(1, numColors, length(idxs))), :));
    %disp(size(times)); disp(size(deformation_amplitudes));
    plot(times, (deformation_amplitudes(idxs, :)), 'LineWidth',2);
    % Only show a subset of 
    set(gca, 'FontSize', 16);
    legend(arrayfun(@(i) string(i), idxs), 'FontSize', 12);
    xlabel('$ t/t_s $', 'Interpreter','Latex', 'FontSize', 20);
    ylabel('$ r/R_o $', 'Interpreter','Latex', 'FontSize', 20);
    sigmaS = default_physical.sigmaS; rhoS = default_physical.rhoS; Ro = default_physical.undisturbed_radius;
    Westar = default_physical.rhoS * default_physical.initial_velocity^2 * ...
        default_physical.undisturbed_radius / default_physical.sigmaS;
    Oh = default_physical.nu * sqrt(default_physical.rhoS / (default_physical.sigmaS ...
        * default_physical.undisturbed_radius));
    Bo = rhoS * default_physical.g * Ro^2 / sigmaS;

    title(sprintf("Deformation Amplitudes with We = %.3g, Oh = %.3g, Bo = %.3g", Westar, Oh, Bo));
    yl = get(gca, 'YLim'); yl = [-max(abs(yl)), max(abs(yl))];
    set(gca, 'YLim', yl); set(gca, 'XLim', 1.05*get(gca, 'XLim'))
  


    saving_figure_pressure = figure('Position', [100, 50, 1000, 300]); % Wider than tall
    hold on;
    %recorded_conditions{ii}.amplitude_defor = 1;
    idxs = 1:(nb_harmonics+1);
    idxs = idxs(max(abs(pressure_amplitudes), [] , 2) < 5); %2.^(1:floor(log2(nb_harmonics)));
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



