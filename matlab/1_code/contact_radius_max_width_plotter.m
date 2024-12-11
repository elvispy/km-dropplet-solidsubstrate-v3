
% Add functions to calculate maximum width
safe_folder = fullfile(fileparts(mfilename('fullpath')), "simulation_code");
addpath(safe_folder, '-begin');


warning('off', 'all');
Ro = 0.0203; % radius in cm
rho = 0.96; %g/cm3
sigma = 20.5; %dyne/m
t_ic = sqrt(rho*Ro^3/sigma); % inertio-capillary time scale
%filename = '../0_data/manual/Low We comparison.xlsx';
filename = '../0_data/manual/Oh comparisons.xlsx';
alldata = load('../2_output/postprocessing.mat', 'data');
alldata = alldata.data;
alldata = alldata(contains(alldata.file_name, 'directComparison') & contains(alldata.parent_folder,'v3') & alldata.number_of_harmonics == 90, :);
sheets = sheetnames(filename);
sheets = sheets(contains(sheets, 'Bounce'));
sheets2 = matlab.lang.makeValidName(sheets);
data = struct();
cmp = "parula";
close all; cmap = parula(100*length(sheets)); ss = size(cmap, 1);
for i = 1:numel(sheets)
    tbl = readtable(filename, 'Sheet', sheets{i}, 'ReadVariableNames', true, 'HeaderLines', 1);
    
    data.(sheets2{i}) = table2struct(tbl, 'ToScalar', true);
    data.(sheets2{i}).We = str2double(regexp(sheets{i}, '(?<=\(We\s*=\s*)\S+(?=\s*\))', 'match', 'once'));

    bnc = data.(sheets2{i});
    % Plotting
    % Plot Contact Time vs Time_s_ (EXPERIMENTAL)
    figure(1); hold on; set(gcf, 'Position', [734 223 ceil(451*16/9) 451]);
    color_vector = repmat(bnc.We, size(bnc.Time_s_));
    idx = ceil(ss*(bnc.We)/4);
    plot(bnc.Time_s_/t_ic, bnc.ContactRadius_mm_/(10*Ro),'LineWidth', (4-bnc.We)/2+1.5, 'DisplayName',"", 'Color', cmap(idx, :));
    %scatter(bnc.Time_s_/t_ic, bnc.ContactRadius_mm_/(10*Ro), 50, cmap(idx, :), 'filled'); %'DisplayName',sprintf("$We=%.2f$", bnc.We));
    
    % Plot Contact Time vs Time_s_ (SIMULATIONS)
    [~, idx2] = min(abs(alldata.weber - bnc.We));  % idx is the index of the closest value
    file = fullfile(pwd, "..", "2_output", alldata.parent_folder(idx2), alldata.file_name(idx2)); 
    values = load(file);
    %alldata = alldata(values.PROBLEM_CONSTANTS.weber == bnc.We, :);
    if  abs(alldata.weber(idx2) - bnc.We) > 1e-3
        disp('couldnt find simulation with value We=', bnc.We);
    else
        length_unit = values.length_unit;
        recorded_conditions = values.recorded_conditions;
        pixel = 0.02*length_unit;%5e-4; % 5e-4 cm %Threshold for experimental contact
        theta_vector = values.PROBLEM_CONSTANTS.theta_vector;
        times_vector_adim = values.recorded_times/t_ic;
        times_simul    = zeros(size(recorded_conditions, 1), 1);
        max_width_adim      = zeros(size(recorded_conditions, 1), 1);
        contact_radius_adim = zeros(size(recorded_conditions, 1), 1);
        drop_top_adim = zeros(size(recorded_conditions, 1), 1);
        drop_top_adim_exp = zeros(size(recorded_conditions, 1), 1);
        drop_bottom_adim = zeros(size(recorded_conditions, 1), 1);
        drop_CM_adim = zeros(size(recorded_conditions, 1), 1);
        for jj = 1:(size(recorded_conditions, 1))
            adim_deformations = recorded_conditions{jj}.deformation_amplitudes/length_unit;
            adim_CM = recorded_conditions{jj}.center_of_mass/length_unit;
            drop_radius = zeta_generator(adim_deformations);
            drop_radius = @(theta) 1 + drop_radius(theta);
            drop_height = @(theta) cos(theta) .* drop_radius(theta) + adim_CM;
            drop_top_adim(jj)    = drop_height(0);
            drop_bottom_adim(jj) = drop_height(pi);
            drop_top_adim_exp(jj)    = max(drop_height(linspace(0, pi/2, 100)));
            drop_CM_adim(jj)     = adim_CM;
                
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
    figure(2); hold on; set(gcf, 'Position', [77 224 ceil(450*17/9) 450]);
    plot(bnc.Time_s_/t_ic, bnc.MaxRadius_mm_/(10*Ro),'Color', cmap(idx, :), 'LineWidth', (4-bnc.We)/2+1.5, 'DisplayName',"");
    %scatter(bnc.Time_s_/t_ic, bnc.MaxRadius_mm_/(10*Ro), 50, cmap(idx, :), 'filled'); %DisplayName',sprintf("$We=%.2f$", bnc.We));
    plot(times_vector_adim, max_width_adim, '--', 'LineWidth', 2, 'DisplayName',"", 'Color', cmap(idx, :));
    
    T = table(times_vector_adim(:)*t_ic, contact_radius_adim(:)*(10*Ro), max_width_adim(:)*(10*Ro), ...
        drop_CM_adim(:)*(10*Ro), drop_bottom_adim(:)*(10*Ro), drop_top_adim(:) * (10*Ro), drop_top_adim_exp(:) * (10*Ro), ...
        'VariableNames', {'Time (s)', 'Contact radius (mm)', 'Max radius (mm)', 'Center of Mass (mm)', 'Bottom (mm)', 'Top (mm)', 'Top (camera view) (mm)'});
        writetable(T, '../2_output/directComparison.xlsx', 'Sheet', sheets{i});
end

figure(1); grid on;
xlabel('$t/t_{ic}$', 'Interpreter','latex');
ylabel('Contact Radius $r/R_o$', 'Interpreter','latex');
%title('Contact Radius vs Time');
%legend('Interpreter','latex');
colormap(cmp);  % You can choose different colormaps (e.g., 'parula', 'jet', 'hot', etc.)
cb = colorbar;  % Show the color scale
caxis([0 3.58]);
ylabel(cb, 'We');
set(gca, 'FontSize', 24);
h1 = plot(NaN, NaN, 'k-', 'MarkerFaceColor', 'k', 'LineWidth', 2);  % Dummy plot for first legend entry
h2 = plot(NaN, NaN, 'k--', 'LineWidth', 2);  % Dummy plot for second legend entry
legend([h1, h2], 'Experiments', 'Simulation');

figure(2); grid on;
xlabel('$t/t_{ic}$', 'Interpreter', ' latex');
ylabel('Maximum Radius ($r/R_o$)', 'Interpreter', 'latex');
%title('Maximum Radius vs Time');
%legend('Interpreter','latex');
colormap(cmp);  % You can choose different colormaps (e.g., 'parula', 'jet', 'hot', etc.)
cb = colorbar;  % Show the color scale
caxis([0 3.58]);
ylabel(cb, 'We');
xlim([0 3]);
set(gca, 'FontSize', 24);
%set(gca,'ColorScale','log')
h1 = plot(NaN, NaN, 'k-', 'MarkerFaceColor', 'k', 'LineWidth', 2);  % Dummy plot for first legend entry
h2 = plot(NaN, NaN, 'k--', 'LineWidth', 2);  % Dummy plot for second legend entry
legend([h1, h2], 'Experiments', 'Simulation')

%filtered_data = rmfield(data, fieldnames(data(cellfun(@(f) ~contains(f, 'Bounce'), fieldnames(data))))');


%x = regexp(str, '(?<=\(We\s*=\s*)\S+(?=\s*\))', 'match', 'once');

function [colors, values] = mapWeToColor(we_query)
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
    we_low = 1e-2;
    we_high = 4;
    c = colorbar;
    % Define log-spaced ticks in the original Oh range
    num_ticks = 5; % Including oh_low and oh_high
    we_ticks = round(logspace(log10(we_low), log10(we_high), num_ticks), 3);
    
    % Map Oh ticks to the normalized range [0, 1] (for the colorbar)
    log_we_low = log10(we_low);
    log_we_high = log10(we_high);
    normalized_ticks = (log10(we_ticks) - log_we_low) / (log_we_high - log_we_low);
    
    % Set the colorbar ticks and labels
    set(c, 'Ticks', normalized_ticks, 'TickLabels', ...
        arrayfun(@(s) sprintf("%.2g", s*(s>we_low)), we_ticks, 'UniformOutput', false), ...
        'FontSize', 16);
    ylabel(c, 'Oh', 'FontSize',18);

    % Validate inputs
    if any(we_query < we_low | we_query > we_high)
        warning('All values in oh_query must be within the range [oh_low, oh_high].');
        we_query(we_query < we_low)  = we_low;
        we_query(we_query > we_high) = we_high;
    end

    % Transform Oh values to logarithmic scale
    log_we_low = log10(we_low);
    log_we_high = log10(we_high);
    log_we_query = log10(we_query);

    % Normalize log-scaled oh_query to [0, 1]
    normalized_oh = (log_we_query - log_we_low) / (log_we_high - log_we_low);

    % Map normalized values to colormap indices
    cmap_size = size(cmap, 1);
    color_indices = round(normalized_oh * (cmap_size - 1)) + 1;

    % Ensure indices are within valid range
    color_indices = max(1, min(cmap_size, color_indices));
    colors = cmap(color_indices, :);
    values = log_we_query;
end
