
% Add functions to calculate maximum width
safe_folder = fullfile(fileparts(mfilename('fullpath')), "simulation_code");
addpath(safe_folder, '-begin');


warning('off', 'all');
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
sampled_values = @(data, bins, per_bin) reshape(randsample(data, bins*per_bin, true, histcounts(data, linspace(min(data), max(data), bins+1))), [], 1);
%expData = sortrows(expData, 'We');
expData = expData(expData.Oh > 0.2 | (mod(1:height(expData), 3)==0)', :);
% Error propagation
g = 9.81;
syms hout vin vout r0 ssigma rrho 
epsilon_error = error_prop((g*(hout-r0) + 0.5 * vout^2)/(0.5* vin^2), [hout, r0, vin vout]);
expData.epsilon_error = epsilon_error(expData.h_t_c__m_,  expData.r_0_m_, ...
    expData.v_0__m_s_, expData.v_t_c__m_s_, expData.m_Pixel, expData.r_0Std, expData.v_0_Stdev, expData.v_t_c_Stdev);
expData.tic = sqrt(expData.m_kg_./expData.V_m3_ .* expData.r_0_m_.^3./expData.sigma_N_m_);
tic_error = error_prop(sqrt(rrho*r0^3/ssigma), [rrho, r0, ssigma]);
expData.tic_error = tic_error(expData.rho_kg_m3_, expData.r_0_m_, expData.sigma_N_m_, ...
    zeros(size(expData.rho_kg_m3_)), expData.r_0Std, zeros(size(expData.sigma_N_m_)));
We_error = error_prop(rrho * vin^2 * r0/ssigma, [rrho, vin, r0, ssigma]);
expData.We_error = We_error(expData.rho_kg_m3_, expData.v_0__m_s_, expData.r_0_m_, expData.sigma_N_m_, ...
    zeros(size(expData.rho_kg_m3_)), ...
    expData.v_0_Stdev,  expData.r_0Std, zeros(size(expData.sigma_N_m_)));
expData.tc_tic = expData.t_c_s_ ./ sqrt(expData.rho_kg_m3_ .* expData.r_0_m_.^3./expData.sigma_N_m_);
alldata = readtable('../2_output/postprocessing.csv');
alldata = alldata(contains(alldata.parent_folder,'v3') & alldata.number_of_harmonics == 90 ...
    & abs(alldata.bond - 0.0189) < 1e-3, :);
alldata = sortrows(alldata, ["weber", "ohnesorge", "bond"]);
Ro = 0.0203;
alldata.tic = sqrt(alldata.bond/981 * Ro);
alldata.ct_adim = 1e-3 * alldata.contact_time_exp_ms./alldata.tic;
% Plotting contact time
figure(1); title('Contact time vs We'); set(gcf, 'Position', [750 176 560 420])
plot_errorbars(expData.We, expData.tc_tic, expData.We_error, expData.tic_error, expData.Oh);
set(gca, 'XScale', 'log', 'FontSize', 16); %, 'YScale', 'log', 'FontSize', 14);
hold on; grid on;
plotVarVsWeByOh(alldata, "ct_adim", ohs, colors);
xticks([0.01, 0.1, 1, 10]);               % Define tick positions
xticklabels({'0.01', '0.1', '1', '10'}); 
xlim([1e-3, 10.5]); ylim([0, 7]);
ylabel('Non-dimensional contac time $t_c/t_{ic}$','Interpreter','latex', 'FontSize',20);
xlabel('We','Interpreter','latex');

% plotting coef of restitution
figure(2); title('Coef res vs We'); set(gcf, 'Position', [50 176 560 420])
[ohs, colors] = plot_errorbars(expData.We, expData.epsilon, expData.We_error, expData.epsilon_error, expData.Oh);
set(gca, 'XScale', 'log', 'FontSize', 16); %, 'YScale', 'log', 'FontSize', 14);
hold on; grid on;
plotVarVsWeByOh(alldata, "coef_rest_exp", ohs, colors);
xticks([0.01, 0.1, 1, 10]);               % Define tick positions
xticklabels({'0.01', '0.1', '1', '10'}); 
xlim([1e-3, 10.5]); ylim([0, 1]);
xlabel('We','Interpreter','latex');
ylabel('Coefficient of restitution ($\varepsilon$)','Interpreter','latex', 'FontSize',20);
return
data = struct();
close all; cmap = jet(10*length(sheets)); ss = size(cmap, 1);
for i = 1:numel(sheets)
    expData = readtable(filename, 'Sheet', sheets{i}, 'ReadVariableNames', true, 'HeaderLines', 2);
    
    data.(sheets{i}) = table2struct(expData, 'ToScalar', true);
    data.(sheets{i}).We = str2double(regexp(sheets{i}, '(?<=\(We\s*=\s*)\S+(?=\s*\))', 'match', 'once'));

    bnc = data.(sheets{i});
    % Plotting
    % Plot Contact Time vs Time_s_ (EXPERIMENTAL)
    figure(1); hold on; set(gcf, 'Position', [734 223 646 451]);
    color_vector = repmat(bnc.We, size(bnc.Time_s_));
    idx = ceil(ss*(bnc.We)/4);
    plot(bnc.Time_s_/t_ic, bnc.ContactRadius_mm_/(10*Ro),'LineWidth', 1, 'DisplayName',"", 'Color', cmap(idx, :));
    scatter(bnc.Time_s_/t_ic, bnc.ContactRadius_mm_/(10*Ro), 50, cmap(idx, :), 'filled'); %'DisplayName',sprintf("$We=%.2f$", bnc.We));
    
    % Plot Contact Time vs Time_s_ (SIMULATIONS)
    [~, idx2] = min(abs(alldata.weber - bnc.We));  % idx is the index of the closest value
    file = fullfile(pwd, "..", "2_output", alldata.parent_folder(idx2), alldata.file_name(idx2)); 
    values = load(file);
    %alldata = alldata(values.PROBLEM_CONSTANTS.weber == bnc.We, :);
    if  abs(alldata.weber - bnc.We) > 1e-3
        disp('couldnt find simulation with value We=', bnc.We);
    else
        length_unit = values.length_unit;
        recorded_conditions = values.recorded_conditions;
        pixel = 5e-4; %Threshold for experimental contact
        theta_vector = values.PROBLEM_CONSTANTS.theta_vector;
        times_vector_adim = values.recorded_times/t_ic;
        times_simul    = zeros(size(recorded_conditions, 1), 1);
        max_width_adim      = zeros(size(recorded_conditions, 1), 1);
        contact_radius_adim = zeros(size(recorded_conditions, 1), 1);
        for jj = 1:(size(recorded_conditions, 1))
            adim_deformations = recorded_conditions{jj}.deformation_amplitudes/length_unit;
            adim_CM = recorded_conditions{jj}.center_of_mass/length_unit;
            drop_radius = zeta_generator(adim_deformations);
            drop_radius = @(theta) 1 + drop_radius(theta);
            drop_height = @(theta) cos(theta) .* drop_radius(theta) + adim_CM;
                
            % Max width calculation & spread time of width
            max_width_adim(jj) = maximum_contact_radius(adim_deformations);
            % if current_width > max_width
            %     max_width = current_width; 
            %     spread_time_width = recorded_times(jj) - touch_time;
            % end
            % EXPERIMENTAL max_contact_radius calculation
            current_contact_points_exp = recorded_conditions{jj}.contact_points;
            if current_contact_points_exp == 0
                current_contact_radius_exp = 0;
            else 
                kk = 1;
                theta2 = theta_vector(current_contact_points_exp+kk);
                while drop_height(theta2) < pixel/length_unit
                    kk = kk + 1;
                    theta2 = theta_vector(current_contact_points_exp + kk);
                end
                theta1 = theta_vector(current_contact_points_exp + kk - 1);
                while abs(theta1-theta2) > pi/1e+3
                    theta_middle = (theta1+theta2)/2;
                    if drop_height(theta_middle) < pixel/length_unit
                        theta1 = theta_middle;
                    else
                        theta2 = theta_middle;
                    end
                end
                current_contact_radius_exp = sin(theta1) ...
                    * drop_radius(theta1);
            end
            contact_radius_adim(jj) = current_contact_radius_exp;
            % if current_contact_radius_exp > max_contact_radius_exp
            %     max_contact_radius_exp = current_contact_radius_exp;
            %     spread_time_exp = recorded_times(jj) - touch_time;
            % end
        end
        plot(times_vector_adim, contact_radius_adim, '--', 'LineWidth', 2, 'DisplayName',"", 'Color', cmap(idx, :));

    end

    % Plot Maximum Radius vs Time_s_
    figure(2); hold on; set(gcf, 'Position', [77 224 1800 450]);
    plot(bnc.Time_s_/t_ic, bnc.MaxRadius_mm_/(10*Ro),'Color', cmap(idx, :), 'LineWidth', 1, 'DisplayName',"");
    scatter(bnc.Time_s_/t_ic, bnc.MaxRadius_mm_/(10*Ro), 50, cmap(idx, :), 'filled'); %DisplayName',sprintf("$We=%.2f$", bnc.We));
    plot(times_vector_adim, max_width_adim, '--', 'LineWidth', 2, 'DisplayName',"", 'Color', cmap(idx, :));
end

figure(1); grid on;
xlabel('$t/t_{ic}$', 'Interpreter','latex');
ylabel('Contact Radius $r/R_o$', 'Interpreter','latex');
title('Contact Radius vs Time');
%legend('Interpreter','latex');
colormap jet;  % You can choose different colormaps (e.g., 'parula', 'jet', 'hot', etc.)
cb = colorbar;  % Show the color scale
caxis([0 3.58]);
ylabel(cb, 'We');
set(gca, 'FontSize', 15);
h1 = plot(NaN, NaN, 'ko-');  % Dummy plot for first legend entry
h2 = plot(NaN, NaN, 'k-');  % Dummy plot for second legend entry
legend([h1, h2], 'Experimental', 'Simulation');

figure(2); grid on;
xlabel('$t/t_{ic}$', 'Interpreter', ' latex');
ylabel('Maximum Radius ($r/R_o$)', 'Interpreter', 'latex');
title('Maximum Radius vs Time');
%legend('Interpreter','latex');
colormap jet;  % You can choose different colormaps (e.g., 'parula', 'jet', 'hot', etc.)
cb = colorbar;  % Show the color scale
caxis([0 3.58]);
ylabel(cb, 'We');
set(gca, 'FontSize', 15);

%filtered_data = rmfield(data, fieldnames(data(cellfun(@(f) ~contains(f, 'Bounce'), fieldnames(data))))');


%x = regexp(str, '(?<=\(We\s*=\s*)\S+(?=\s*\))', 'match', 'once');

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

    % Create scatter plot
    scatter(x, y, 50, z, 'filled'); % Size 100, color by z
    colormap("winter"); % Choose a colormap
    c = colorbar; % Add a colorbar
    ylabel(c, 'Oh');
    hold on;

    % Get the colormap and normalize z to map colors correctly
    cmap = colormap;
    z_norm = (z - min(z)) / (max(z) - min(z)); % Normalize z to [0, 1]
    colors = interp1(linspace(0, 1, size(cmap, 1)), cmap, z_norm);
    [zs, iid] = sort(z_norm);
    unique_colors = interp1(linspace(0, 1, size(cmap, 1)), cmap, zs);
    
    % Add error bars with the same color as scatter points
    for i = 1:length(x)
        errorbar(x(i), y(i), ey(i), ey(i), ex(i), ex(i), 'LineWidth', 1, ...
            'Color', colors(i, :))
        % Horizontal error bar
        %line([x(i)-ex(i), x(i)+ex(i)], [y(i), y(i)], ...
        %     'Color', colors(i, :), 'LineWidth', 1);
        % Vertical error bar
        %line([x(i), x(i)], [y(i)-ey(i), y(i)+ey(i)], ...
        %     'Color', colors(i, :), 'LineWidth', 1);
    end

    hold off;
    unique_colors = unique_colors(iid, :);
end

function plotVarVsWeByOh(table, var, ohs, colors)
    % Get unique Oh values from the table
    unique_ohs = unique(table.ohnesorge);

    % Loop through each unique Oh value
    for i = 1:length(unique_ohs)
        target_oh = unique_ohs(i);
        
        % Find the closest Oh index
        [~, idx] = min(abs(ohs - target_oh));  
        selected_color = colors(idx, :);  % Get corresponding color
        
        % Filter table for current Oh value
        filtered_table = table(abs(table.ohnesorge - target_oh) < 1e-6, :);  % Adjust tolerance if needed
        
        % Plot epsilon vs We
        plot(filtered_table.weber, filtered_table.(var), 'Color', selected_color, ...
            'LineWidth', 3, 'LineStyle','--');
        hold on;
    end

end
