
% Add functions to calculate maximum width
safe_folder = fullfile(fileparts(mfilename('fullpath')), "simulation_code");
addpath(safe_folder, '-begin');


%warning('off', 'all');
close all;
Ro = 0.0203; % radius in cm
rho = 0.96; %g/cm3
sigma = 20.5; %dyne/m
t_ic = sqrt(rho*Ro^3/sigma); % inertio-capillary time scale
filename = '../0_data/manual/CleanDataAll.xlsx';

sheets = sheetnames(filename);
%sheets = sheets(contains(sheets, 'Bounce'));
%sheets = matlab.lang.makeValidName(sheets);
expData = readtable(filename, 'Sheet', sheets, 'ReadVariableNames', true, 'HeaderLines', 2);
%sampled_values = @(data, bins, per_bin) reshape(randsample(data, bins*per_bin, true, histcounts(data, linspace(min(data), max(data), bins+1))), [], 1);
%expData = sortrows(expData, 'We');
%expData = expData(expData.Oh > 0.2 | (mod(1:height(expData), 1)==0)', :);
% Error propagation
g = 9.81;
syms hout vin vout r0 ssigma rrho tc
epsilon_error = error_prop(sqrt((g*(hout-r0) + 0.5 * vout^2)/(0.5* vin^2)), [hout, r0, vin vout]);
expData.epsilon_error = epsilon_error(expData.h_t_c__m_,  expData.r_0_m_, ...
    expData.v_0__m_s_, expData.v_t_c__m_s_, expData.m_Pixel, expData.r_0Std, expData.v_0_Stdev, expData.v_t_c_Stdev);
expData.tic = sqrt(expData.m_kg_./expData.V_m3_ .* expData.r_0_m_.^3./expData.sigma_N_m_);
tic = sqrt(expData.rho_kg_m3_ .* expData.r_0_m_.^3./expData.sigma_N_m_);
tc_tic_error = error_prop(tc/sqrt(rrho*r0^3/ssigma), [tc, rrho, r0, ssigma]);
expData.tc_tic_error = tc_tic_error(expData.t_c_s_, expData.rho_kg_m3_, expData.r_0_m_, expData.sigma_N_m_, ...
    expData.deltaT_s_, zeros(size(expData.rho_kg_m3_)), expData.r_0Std, zeros(size(expData.sigma_N_m_)));
We_error = error_prop(rrho * vin^2 * r0/ssigma, [rrho, vin, r0, ssigma]);
expData.We_error = We_error(expData.rho_kg_m3_, expData.v_0__m_s_, expData.r_0_m_, expData.sigma_N_m_, ...
    zeros(size(expData.rho_kg_m3_)), ...
    expData.v_0_Stdev,  expData.r_0Std, zeros(size(expData.sigma_N_m_)));
expData.tc_tic = expData.t_c_s_ ./ tic;

%% Reading linear model data
alldata = readtable('../2_output/postprocessing.csv');
alldata = alldata(contains(alldata.parent_folder,'v3') & alldata.number_of_harmonics == 90 ...
    & abs(alldata.bond - 0.0189) < 1e-3, :);
alldata = sortrows(alldata, ["weber", "ohnesorge", "bond"]);
Ro = 0.0203;
alldata.tic = sqrt(alldata.bond/981 * Ro);
alldata.ct_adim = 1e-3 * alldata.contact_time_exp_ms./alldata.tic;

%% Reading DNS
filename = '../0_data/manual/Drop Parameter spreadsheet.xlsx';
sheets = sheetnames(filename);
sheets = sheets(contains(sheets, 'Level 12'));
%sheets = matlab.lang.makeValidName(sheets);
DNSData = readtable(filename, 'Sheet', sheets, 'ReadVariableNames', true, 'HeaderLines', 1);
DNSData.tic = sqrt(DNSData.rho_kg_m3_ .* DNSData.r_0_m_.^3 ./ DNSData.sigma_N_m_);
DNSData.tc_tic = DNSData.t_c_s_ ./ DNSData.tic;

%% Plotting

% plotting coef of restitution
f2 = figure(2); set(gcf, 'Position', [5 176 ceil(420*18/9) 420])
colormap("parula"); % Choose a colormap
plot_errorbars(expData.We, expData.epsilon, expData.We_error, expData.epsilon_error, expData.Oh);
hold on; grid on;
plotVarVsWeByOh(alldata, "coef_rest_exp");
scatterDNS(DNSData, "epsilon");
%custom_ticks = 0.2:0.2:0.8; % Specify desired tick positions
%set(c, 'Ticks', custom_ticks, 'TickLabels', arrayfun(@num2str, custom_ticks, 'UniformOutput', false));
set(gca, 'XScale', 'log', 'FontSize', 16); %, 'YScale', 'log', 'FontSize', 14);
xticks([0.01, 0.1, 1, 10]);               % Define tick positions
xticklabels({'0.01', '0.1', '1', '10'}); 
xlim([1e-3, 12]); ylim([0, 1]);
xlabel('$We = \rho V_0^2 R_0 / \sigma$','Interpreter','latex');
ylabel('Coef. Restitution ($\varepsilon = \sqrt{E_{out}/E_{in}}$)','Interpreter','latex', 'FontSize',20);
h1 = scatter(nan, nan, 250, "^", "filled", ...
            'MarkerFaceColor', 'black', 'MarkerEdgeColor','k', 'LineWidth',1.5);
h2 = plot(NaN, NaN, 'k-', 'LineWidth', 2);  % Dummy plot for second legend entry
h3 = plot(NaN, NaN, 'ko', 'LineWidth', 2);  % Dummy plot for second legend entry
legend([h1, h2, h3], 'DNS', 'Kinematic Match', 'Experiments', 'Location','southwest', 'FontSize', 18);
saveas(f2, fullfile(safe_folder, "..", "..", "2_output", "Figures", "epsilonvsWeExp+DNS.png"));

% Plotting contact time
f1 = figure(1); set(gcf, 'Position', [950 176 ceil(420*13/9) 420])
colormap("parula"); % Choose a colormap
plot_errorbars(expData.We, expData.tc_tic, expData.We_error, expData.tc_tic_error, expData.Oh);
set(gca, 'XScale', 'log', 'FontSize', 16); %, 'YScale', 'log', 'FontSize', 14);
hold on; grid on;
plotVarVsWeByOh(alldata, "ct_adim");
scatterDNS(DNSData, "tc_tic");
xticks([0.01, 0.1, 1, 10]);               % Define tick positions
xticklabels({'0.01', '0.1', '1', '10'}); 
xlim([1e-3, 12]); ylim([0, 7]);
ylabel('$t_c/t_{ic}$','Interpreter','latex', 'FontSize',20);
xlabel('$We = \rho V_0^2 R_0 / \sigma$','Interpreter','latex');
h1 = scatter(nan, nan, 250, "^", "filled", ...
            'MarkerFaceColor', 'black', 'MarkerEdgeColor','k', 'LineWidth',1.5);
h2 = plot(NaN, NaN, 'k-', 'LineWidth', 2);  % Dummy plot for second legend entry
h3 = plot(NaN, NaN, 'ko', 'LineWidth', 2);  % Dummy plot for second legend entry
legend([h1, h2, h3], 'DNS', 'Kinematic Match', 'Experiments', 'Location','southwest', 'FontSize', 18);
saveas(f1, fullfile(safe_folder, "..", "..", "2_output", "Figures", "tcvsWeExp+DNS.png"));


%% Aux functions
function sigma_f_func = error_prop(expr, variables)
    % Create symbolic variables for uncertainties
    uncertainties = sym("s" + string(variables));
    
    % Compute the uncertainty propagation formula
    terms = arrayfun(@(var, svar) diff(expr, var)^2 * svar^2, ...
                     variables, uncertainties, 'UniformOutput', false);
    sigma_f = sqrt(sum(cat(1, terms{:}))); % Concatenate cell contents into an array
    % Combine variables and their uncertainties
    all_vars = [variables, uncertainties];
    
    % Generate a MATLAB function using matlabFunction
    sigma_f_func = matlabFunction(sigma_f, 'Vars', all_vars);
end



function [zs, unique_colors] = plot_errorbars(x, y, ex, ey, z)
    % Function to plot scatter points with error bars and color mapping
    %
    % Inputs:
    %   x  - Vector of x-coordinates
    %   y  - Vector of y-coordinates
    %   ex - Vector of x-error bar values
    %   ey - Vector of y-error bar values
    %   z  - Vector of values to color the scatter points and error bars
    %
    % Example:
    %   plot_errorbars([1, 2, 3], [2, 3, 1], [0.1, 0.2, 0.1], [0.2, 0.1, 0.2], [10, 20, 30]);

    hold on;
    zs = nan; unique_colors = nan;
    % Get the colormap and normalize z to map colors correctly
    % z_trans = log10(z+1e-2 * (z <=0));
    % z_norm = (z_trans - min(z_trans)) / (max(z_trans) - min(z_trans)); % Normalize z to [0, 1]
    % colors = interp1(linspace(0, 1, size(cmap, 1)), cmap, z_norm);
    % [zs, iid] = sort(z_norm);
    % unique_colors = interp1(linspace(0, 1, size(cmap, 1)), cmap, zs);
    % unique_colors = unique_colors(iid, :);
    % 
    % Create scatter plot
    colors = mapOhToColor(z);
    scatter(x, y, 20, colors, 'filled'); % Size 100, color by z
    % Add error bars with the same color as scatter points
    %return
    for i = 1:length(x)
        if 2*ey(i) >0.6; continue; end 
        errorbar(x(i), y(i), ey(i), ey(i), ex(i), ex(i), 'LineWidth', 1, ...
            'Color', colors(i, :))
    end

    hold off;
    
end

function plotVarVsWeByOh(table, var)
    % Get unique Oh values from the table
    unique_ohs = unique(table.ohnesorge);

    % Loop through each unique Oh value
    for i = 1:length(unique_ohs)
        target_oh = unique_ohs(i);
        target_color = mapOhToColor(target_oh);
        % Find the closest Oh index
        %[~, idx] = min(abs(ohs - target_oh));  
        %selected_color = colors(idx, :);  % Get corresponding color
        
        % Filter table for current Oh value
        filtered_table = table(table.ohnesorge == target_oh, :);  % Adjust tolerance if needed
        
        % Plot epsilon vs We
        plot(filtered_table.weber, filtered_table.(var), 'Color', target_color, ...
            'LineWidth', 5, 'LineStyle','-');
        hold on;
    end

end

function scatterDNS(table, var)
    % Get unique Oh values from the table
    unique_ohs = unique(table.Oh);

    % Loop through each unique Oh value
    for i = 1:length(unique_ohs)
        target_oh = unique_ohs(i);
        target_color = mapOhToColor(target_oh);
        % Find the closest Oh index
        %[~, idx] = min(abs(ohs - target_oh));  
        %selected_color = colors(idx, :);  % Get corresponding color
        
        % Filter table for current Oh value
        filtered_table = table(abs(table.Oh - target_oh) < 1e-4, :);  % Adjust tolerance if needed
        
        % Plot epsilon vs We
        scatter(filtered_table.We, filtered_table.(var), 250, "^", "filled", ...
            'MarkerFaceColor', target_color, 'MarkerEdgeColor','k', 'LineWidth',1.5);
        hold on;
    end

end

function [colors, values] = mapOhToColor(oh_query)
    % Maps a vector of Ohnesorge numbers (oh_query) to indices in the colormap (cmap)
    % on a logarithmic scale between oh_low and oh_high.
    %
    % Inputs:
    %   oh_low   - Minimum Ohnesorge value (lower bound)
    %   oh_high  - Maximum Ohnesorge value (upper bound)
    %   oh_query - Vector of Ohnesorge values to query
    %   cmap     - Colormap matrix (Nx3)
    %
    % Output:
    %   color_indices - Array of indices corresponding to oh_query in cmap
    
    cmap = colormap;
    oh_low = 1e-3;
    oh_high = 0.8;
    c = colorbar;
    % Define log-spaced ticks in the original Oh range
    num_ticks = 5; % Including oh_low and oh_high
    oh_ticks = round(logspace(log10(oh_low), log10(oh_high), num_ticks), 3);
    
    % Map Oh ticks to the normalized range [0, 1] (for the colorbar)
    log_oh_low = log10(oh_low);
    log_oh_high = log10(oh_high);
    normalized_ticks = (log10(oh_ticks) - log_oh_low) / (log_oh_high - log_oh_low);
    
    % Set the colorbar ticks and labels
    set(c, 'Ticks', normalized_ticks, 'TickLabels', ...
        arrayfun(@(s) sprintf("%.2g", s*(s>oh_low)), oh_ticks, 'UniformOutput', false), ...
        'FontSize', 16);
    ylabel(c, 'Oh', 'FontSize',18);

    % Validate inputs
    if any(oh_query < oh_low | oh_query > oh_high)
        warning('All values in oh_query must be within the range [oh_low, oh_high].');
        oh_query(oh_query < oh_low)  = oh_low;
        oh_query(oh_query > oh_high) = oh_high;
    end

    % Transform Oh values to logarithmic scale
    log_oh_low = log10(oh_low);
    log_oh_high = log10(oh_high);
    log_oh_query = log10(oh_query);

    % Normalize log-scaled oh_query to [0, 1]
    normalized_oh = (log_oh_query - log_oh_low) / (log_oh_high - log_oh_low);

    % Map normalized values to colormap indices
    cmap_size = size(cmap, 1);
    color_indices = round(normalized_oh * (cmap_size - 1)) + 1;

    % Ensure indices are within valid range
    color_indices = max(1, min(cmap_size, color_indices));
    colors = cmap(color_indices, :);
    values = log_oh_query;
end
