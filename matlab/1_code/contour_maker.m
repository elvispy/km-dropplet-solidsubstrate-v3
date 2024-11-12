

% Load the table
data = readtable('../2_output/postprocessing.csv');  % Replace with the actual file name

% Filter rows
data = data(contains(data.file_name, 'Sweep'), :);


% Extract relevant columns and remove rows with NaNs in the output
weber = data.weber;
ohnesorge = data.ohnesorge;
output_data = data.coef_rest_exp;
bond = data.bond;

% Remove rows where contact_time_ms has NaN values
valid_indices = (~isnan(output_data) & ohnesorge >= 0 & bond == 0);
weber = weber(valid_indices); 
ohnesorge = ohnesorge(valid_indices);
ohnesorge(ohnesorge == 0) = 1e-9;
output_data = output_data(valid_indices);

% Create a fine grid in log-space for interpolation
log_weber = logspace(log10(min(weber)), log10(max(weber)), 100);
lin_ohnesorge = logspace(log10(min(ohnesorge)), log10(max(ohnesorge)), 100);
[WeberGrid, OhnesorgeGrid] = meshgrid(log_weber, lin_ohnesorge);

% Interpolate the contact time on this grid
ContactTimeGrid = griddata(weber, ohnesorge, output_data, WeberGrid, OhnesorgeGrid, 'cubic');

% Plot contour with log scales
contourf(WeberGrid, OhnesorgeGrid, ContactTimeGrid, 20);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 14);
xlabel('Weber', 'FontSize',16);
ylabel('Ohnesorge', 'FontSize',16);
yticks = [0.02:0.04:0.1, 0.2:0.4:1]; 
ylim([0.015, max(ohnesorge)])
set(gca, 'YTick', yticks, 'YTickLabel', yticks);  % Apply linear tick labels
% Customize x-axis tick appearance
set(gca, 'TickDir', 'in', 'TickLength', [0.03, 0.01]);  % Increase tick length on x-axis

colorbar;
title('Contour plot $\varepsilon_{exp}$', 'Interpreter','latex', 'FontSize',24);