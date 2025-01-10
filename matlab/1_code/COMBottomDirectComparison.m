
% Add functions to calculate maximum width
safe_folder = fullfile(fileparts(mfilename('fullpath')), "simulation_code");
addpath(safe_folder, '-begin');


warning('off', 'all');
%Ro = 0.0203; % radius in cm
%rho = 0.87; %g/cm3
%sigma = 18.7; %dyne/m
%t_ic = sqrt(rho*Ro^3/sigma); % inertio-capillary time scale
%filename = '../0_data/manual/Low We comparison.xlsx';
%filename = '../0_data/manual/Oh comparisons.xlsx';
alldata = load('../2_output/postprocessing.mat', 'data');
alldata = alldata.data;
alldata = alldata(contains(alldata.file_name, 'directComparison') & contains(alldata.parent_folder,'v3') & alldata.number_of_harmonics == 90, :);
%sheets = sheetnames(filename);
%sheets = sheets(contains(sheets, 'Bounce'));
%sheets2 = matlab.lang.makeValidName(sheets);
data = struct();
%cmp = "parula";
%close all; cmap = parula(100*length(sheets)); ss = size(cmap, 1);
root_folder = fileparts(fileparts(mfilename('fullpath')));
files_folder = dir(fullfile(root_folder, "2_output", "**/*.mat"));

% Values of nondimensional 
Ohs = [0.0303, 0.0303, 0.0298]';
Wes = [0.0231, 0.253, 2.269]';
Bos = [0.019, 0.019, 0.020153]';

for ii = 1:length(Ohs)
    %tbl = readtable(filename, 'Sheet', sheets{ii}, 'ReadVariableNames', true, 'HeaderLines', 1);
    
    %load(alldata.)
    %data(ii) = table2struct(tbl, 'ToScalar', true);
    data(ii).We = Wes(ii);
    data(ii).Oh = Ohs(ii);
    data(ii).Bo = Bos(ii);

    bnc = data(ii);
    % Plotting
    % Plot Contact Time vs Time_s_ (EXPERIMENTAL)
    %figure(1); hold on; set(gcf, 'Position', [734 223 ceil(451*16/9) 451]);
    %color_vector = repmat(bnc.We, size(bnc.Time_s_));
    %idx = ceil(ss*(bnc.We)/4);
    %plot(bnc.Time_s_/t_ic, bnc.ContactRadius_mm_/(10*Ro),'LineWidth', (4-bnc.We)/2+1.5, 'DisplayName',"", 'Color', cmap(idx, :));
    %scatter(bnc.Time_s_/t_ic, bnc.ContactRadius_mm_/(10*Ro), 50, cmap(idx, :), 'filled'); %'DisplayName',sprintf("$We=%.2f$", bnc.We));
    
    % Plot Contact Time vs Time_s_ (SIMULATIONS)
    alldata2 = alldata(abs(alldata.weber - bnc.We)/bnc.We < 5e-3 &  ...
        abs(alldata.ohnesorge - bnc.Oh)/bnc.Oh < 5e-3 &  ...
        abs(alldata.bond - bnc.Bo)/bnc.Bo < 5e-3, :);
    %[~, idx2] = min(abs(alldata2.weber - bnc.We));  % idx is the index of the closest value
    
    %alldata = alldata(values.PROBLEM_CONSTANTS.weber == bnc.We, :);
    if  height(alldata2) ~= 1
        disp('couldnt find simulation with value We=', bnc.We);
    else
        file = fullfile(pwd, "..", "2_output", alldata2.parent_folder(1), alldata2.file_name(1)); 
        values = load(file);
        length_unit = values.length_unit;
        recorded_conditions = values.recorded_conditions;
        Ro = values.default_physical.undisturbed_radius;
        tic = sqrt(values.default_physical.rhoS*values.default_physical.undisturbed_radius.^3 ...
            /values.default_physical.sigmaS); % inertio-capillary time scale
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
        M = 200;
        dropshape = zeros(size(recorded_conditions, 1), M);
        for jj = 1:(size(recorded_conditions, 1))
            adim_deformations = recorded_conditions{jj}.deformation_amplitudes/length_unit;
            adim_CM = recorded_conditions{jj}.center_of_mass/length_unit;
            drop_radius = zeta_generator(adim_deformations);
            drop_radius = @(theta) 1 + drop_radius(theta);
            dropshape(jj, :) = drop_radius(linspace(0, pi, M))*length_unit;
            drop_height = @(theta) cos(theta) .* drop_radius(theta) + adim_CM;
            drop_top_adim(jj)    = drop_height(0);
            drop_bottom_adim(jj) = drop_height(pi);
            drop_top_adim_exp(jj)= max(drop_height(linspace(0, pi/2, 100)));
            drop_CM_adim(jj)     = adim_CM;
                
            % Max width calculation & spread time of width
            max_width_adim(jj) = maximum_contact_radius(adim_deformations);
            % if current_width > max_width
            %     max_width = current_width; 
            %     spread_time_width = recorded_times(jj) - touch_time;
            % end
            % EXPERIMENTAL max_contact_radius calculation
            current_contact_points_exp = recorded_conditions{jj}.contact_points;
            if current_contact_points_exp == 0 && false
                current_contact_radius_exp = 0;
            else 
                kk = 1;
                theta2 = theta_vector(current_contact_points_exp+kk);
                while drop_height(theta2) < pixel/length_unit
                    kk = kk + 1;
                    theta2 = theta_vector(current_contact_points_exp + kk);
                end
                theta1 = theta_vector(max(current_contact_points_exp + kk - 1, 1));
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
            
        end
        %plot(times_vector_adim, contact_radius_adim, '--', 'LineWidth', 2, 'DisplayName',"", 'Color', cmap(idx, :));
       
        V0 = abs(values.default_physical.initial_velocity);
        g = values.default_physical.g;
        t0 = (-V0 + sqrt(V0^2 - 2*g*0.02*Ro))/g; % Calculating experimental start of contact to substract
        X = @(t) -t * V0 - g*t.^2 + Ro; 
        rc = @(t) 10* sqrt(Ro^2 - (X(t) - 0.02*Ro).^2); % 10 because we go from cm to mm
    end
    tts = linspace(0, -0.8*t0, 5)';
    times = [tts ;times_vector_adim(:)*t_ic - t0];
    c_radii = [rc(tts); contact_radius_adim(:)*(10*Ro)];
    max_width = [10*Ro*ones(5, 1); max_width_adim(:)*(10*Ro)];
    CM = [10*X(tts); drop_CM_adim(:)*(10*Ro)];
    bottom =  [10*(X(tts)-Ro);drop_bottom_adim(:)*(10*Ro)];
    top    =  [10*(X(tts)+Ro);drop_top_adim(:)*(10*Ro)];
    top_exp = [10*(X(tts)+Ro);drop_top_adim_exp(:) * (10*Ro)];

    T = table(times, c_radii, max_width, CM, bottom, top, top_exp, ...
        'VariableNames', {'Time (s)', 'Contact radius (mm)', 'Max radius (mm)', 'Center of Mass (mm)', 'Bottom (mm)', 'Top (mm)', 'Top (camera view) (mm)'});
        writetable(T, '../2_output/directComparisonNew.xlsx', 'Sheet', sprintf("We=%.5g", bnc.We));
    writematrix([[nan, linspace(0, pi, M)]; ...
        [times, [repmat(dropshape(1, :), 5, 1); dropshape]]], ...
        '../2_output/directComparisonShape.xlsx', 'Sheet', sprintf("We=%.5g", bnc.We));
end
