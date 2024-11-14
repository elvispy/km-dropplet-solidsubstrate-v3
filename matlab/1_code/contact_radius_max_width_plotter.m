
% Add functions to calculate maximum width
safe_folder = fullfile(fileparts(mfilename('fullpath')), "simulation_code");
addpath(safe_folder, '-begin');


warning('off', 'all');
Ro = 0.0203; % radius in cm
rho = 0.96; %g/cm3
sigma = 20.5; %dyne/m
t_ic = sqrt(rho*Ro^3/sigma); % inertio-capillary time scale
filename = '../0_data/manual/Low We comparison.xlsx';
alldata = load('../2_output/postprocessing.mat', 'data');
alldata = alldata.data;
alldata = alldata(abs(alldata.ohnesorge - 0.03) < 0.05, :);
sheets = sheetnames(filename);
sheets = sheets(contains(sheets, 'Bounce'));
sheets2 = matlab.lang.makeValidName(sheets);
data = struct();
close all; cmap = jet(10*length(sheets)); ss = size(cmap, 1);
for i = 1:numel(sheets)
    tbl = readtable(filename, 'Sheet', sheets{i}, 'ReadVariableNames', true, 'HeaderLines', 1);
    
    data.(sheets2{i}) = table2struct(tbl, 'ToScalar', true);
    data.(sheets2{i}).We = str2double(regexp(sheets{i}, '(?<=\(We\s*=\s*)\S+(?=\s*\))', 'match', 'once'));

    bnc = data.(sheets2{i});
    % Plotting
    % Plot Contact Time vs Time_s_ (EXPERIMENTAL)
    figure(1); hold on; set(gcf, 'Position', [734 223 646 451]);
    color_vector = repmat(bnc.We, size(bnc.Time_s_));
    idx = ceil(ss*(bnc.We)/4);
    plot(bnc.Time_s_/t_ic, bnc.ContactRadius_mm_/(10*Ro),'LineWidth', 1, 'DisplayName',"", 'Color', cmap(idx, :));
    scatter(bnc.Time_s_/t_ic, bnc.ContactRadius_mm_/(10*Ro), 50, cmap(idx, :), 'filled'); %'DisplayName',sprintf("$We=%.2f$", bnc.We));
    
    % Plot Contact Time vs Time_s_ (SIMULATIONS)
    [~, idx2] = min(abs(alldata.weber - bnc.We));  % idx is the index of the closest value
    file = fullfile(pwd, "..", "2_output", alldata.parent_folder(idx), alldata.file_name(idx2)); 
    values = load(file);
    %alldata = alldata(values.PROBLEM_CONSTANTS.weber == bnc.We, :);
    if values.PROBLEM_CONSTANTS.weber ~= bnc.We
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
        for jj = 1:(size(recorded_conditions, 1)-1)
            adim_deformations = recorded_conditions{jj}.deformation_amplitudes/length_unit;
            adim_CM = recorded_conditions{jj}.center_of_mass/length_unit;
            drop_radius = zeta_generator(adim_deformations);
            drop_radius = @(theta) 1 + drop_radius(theta);
            drop_height = @(theta) cos(theta) .* drop_radius(theta) + adim_CM;

            
            % % contact_time and coef_res
            % if isnan(touch_time)
            %     if recorded_conditions{jj}.contact_points == 0 && ...
            %             recorded_conditions{jj+1}.contact_points > 0
            %         touch_time = recorded_times(jj);
            %         Vin = recorded_conditions{jj}.center_of_mass_velocity;
            %         CM_in = recorded_conditions{jj}.center_of_mass;
            %         Ein = 1/2 * Vin^2;
            %     end
            % end
            
            % if isnan(liftoff_time)
            %     if recorded_conditions{jj}.contact_points > 0  && ...
            %             recorded_conditions{jj+1}.contact_points == 0
            %         liftoff_time = recorded_times(jj);
            %         Vout = recorded_conditions{jj}.center_of_mass_velocity;
            %         Eout = 1/2*Vout^2 + (recorded_conditions{jj}.center_of_mass ...
            %             - CM_in)*g;
            %     end
            % end
            % if ~isnan(touch_time) && ~isnan(liftoff_time)
            %     contact_time = liftoff_time - touch_time;
            %     coef_restitution = sqrt(abs(Eout/Ein));
            % end

            % Min_height (north pole)
            % current_height = drop_height(0);
            % next_height = drop_height(pi/100);
            % if current_height < north_pole_min_height; north_pole_min_height = current_height; end
            % if next_height < current_height %&& current_height < north_pole_exp_min_height
            %     %north_pole_exp_min_height = current_height;
            %     if current_height < north_pole_exp_min_height; north_pole_exp_min_height = current_height; end
            % else
            %     currang = pi/200; dth = pi/200;
            %     next_next_height = drop_height(currang + dth);
            %     while next_next_height >= next_height && currang < pi/3
            %         currang = currang + dth;
            %         next_height = next_next_height;
            %         next_next_height = drop_height(currang + dth);
            %     end
            %     if next_height < north_pole_exp_min_height; north_pole_exp_min_height = next_height; end
            % end

            % % Experimental contact_radius
            % if isnan(liftoff_time_exp)
            %     % If contact ended numerically but simulation ended, record lift off time anyways
            %     if drop_height(pi) > pixel/length_unit %||((size(recorded_conditions, 1)-1 == jj && ~isnan(liftoff_time)))
            %         liftoff_time_exp = recorded_times(jj);
            %         Vout_exp = recorded_conditions{jj}.center_of_mass_velocity;
            %         Eout_exp = 1/2*Vout_exp^2 + (recorded_conditions{jj}.center_of_mass ...
            %             - CM_in)*g;
            %         % if (size(recorded_conditions, 1)-1 == jj && ~isnan(liftoff_time)) 
            %         %     warning("Simulation ended abruptly but numerical contact has ended. We recorded the experimental end of contact anyways for simul \n %s", files_folder(ii).name);
            %         % end
            % 
            %     end
            % end
                
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
    figure(2); hold on; set(gcf, 'Position', [77 224 600 452]);
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
clim([0 3.58]);
ylabel(cb, 'We');
set(gca, 'FontSize', 15);
h1 = plot(NaN, NaN, 'ko-');  % Dummy plot for first legend entry
h2 = plot(NaN, NaN, 'k-');  % Dummy plot for second legend entry
legend([h1, h2], 'Experimental', 'Simulation');

figure(2); grid on;
xlabel('Time (s)');
ylabel('Maximum Radius (mm)');
title('Maximum Radius vs Time');
%legend('Interpreter','latex');
colormap jet;  % You can choose different colormaps (e.g., 'parula', 'jet', 'hot', etc.)
cb = colorbar;  % Show the color scale
clim([0 3.58]);
ylabel(cb, 'We');
set(gca, 'FontSize', 15);

%filtered_data = rmfield(data, fieldnames(data(cellfun(@(f) ~contains(f, 'Bounce'), fieldnames(data))))');


%x = regexp(str, '(?<=\(We\s*=\s*)\S+(?=\s*\))', 'match', 'once');
