function We = criticalWeber(Oh, Bo)
% Oh = 0.0151883524683513;
% Bo = 0.0189312135804878;
rho = 0.96; sigma = 20.5; g = 981;
Ro = sqrt(Bo * sigma / (rho *g));
nu = Oh .* sqrt(sigma*Ro/rho);
harmonics_qtt = 90; version = 3;

numerical_parameters = struct("harmonics_qtt", 90,...
            "simulation_time", inf, "version", 3);
options = struct('version', 3, 'save_results', false, 'optimize_for_bounce', true);
physical_parameters = struct("undisturbed_radius", Ro, ...
            "initial_height", nan, "initial_velocity", nan, ...
            "initial_amplitudes", zeros(1, 90), ...
            "pressure_amplitudes", zeros(1, 90+1), "initial_contact_points", 0, ...
            "rhoS", rho, "sigmaS", sigma, "nu", nu, "g", g);

[guesses, matrix] = load_previous_guesses(harmonics_qtt, version); 
idx = find(and(matrix.Oh == Oh, matrix.Bo == Bo));
if ~isempty(idx)
    flag = false;
    We = matrix{idx(1), 'We'};
else
    [We, flag] = next_guess(guesses, Oh, Bo);
end

while flag == true
    physical_parameters.initial_velocity = - sqrt(We .* sigma./(rho * Ro));
    [recorded_conditions, ~, ~] = ...
        solve_motion_v2(physical_parameters, numerical_parameters, options);
    bounce = is_there_bounce(recorded_conditions, Ro);
    guesses(end+1, :) = {bounce, We, Oh, Bo, ''};
    [We, flag] = next_guess(guesses, Oh, Bo);
end

matrix(end+1, :) = {We, Oh, Bo};
writetable(matrix, '../2_output/criticalWeEvals.csv');
writetable(guesses, sprintf('../2_output/bouncedataharm_qtt=%dversion%d.csv', harmonics_qtt, version));

% Load Python3 in MACOS based on https://www.mathworks.com/matlabcentral/answers/359408-trouble-with-a-command-in-matlab-s-system
if ~ispc && system('python3 --version') ~= 0; setenv('PATH', [getenv('PATH') ':/usr/local/bin/']); end

system('python3 sending_email.py'); % Sending email to notify that's finished
end % end main function

function [data, matrix] = load_previous_guesses(harmonics_qtt, version)
    try
        matrix = readtable('../2_output/criticalWeEvals.csv');
    catch
        %warning('couldnt find previous critical weber evaluations');
        matrix = array2table(zeros(0,3),'VariableNames', {'We', 'Oh', 'Bo'});
    end

    % FInding all .mat files that correspond to simulations
    root_folder = fileparts(fileparts(mfilename('fullpath')));
    files_folder = dir(fullfile(root_folder, "2_output", "**/*.mat"));

    data = array2table(nan(length(files_folder),5),...
        'VariableNames', {'bounce', 'We', 'Oh', 'Bo', 'file_name'});
    data.Properties.VariableTypes(end) = "string";
    %'VariableTypes', ["double", "double", "double", "double", "string"]);
    try
        data2 = readtable(sprintf('../2_output/bouncedataharm_qtt=%dversion%d.csv', harmonics_qtt, version));
        data(1:height(data2), :) = data2;

    catch
        nan;
    end
    for ii = 1:length(files_folder)
        
        if contains(files_folder(ii).name, "postprocessing") || ...
                contains(lower(files_folder(ii).name), "error")  || ...
                (ismember(files_folder(ii).name, data.file_name) && ~isempty(char(files_folder(ii).name)))
            continue; 
        end
        %if ~isnan(data{ii, "coef_rest_exp"}); continue; end
        lastwarn('', ''); clear recorded_conditions recorded_times default_physical length_unit theta_vector
        load(fullfile(files_folder(ii).folder, files_folder(ii).name), ...
            "recorded_conditions", "default_physical", ...
            'default_numerical', "length_unit");
        if contains(lastwarn, 'not found'); error(lastwarn);  end

        % Define parameters of bounce 
        Ro = default_physical.undisturbed_radius;
        We = default_physical.rhoS * default_physical.initial_velocity^2 * ...
            Ro / default_physical.sigmaS;
        Oh = default_physical.nu * sqrt(default_physical.rhoS / (default_physical.sigmaS ...
            * default_physical.undisturbed_radius));
        Bo = default_physical.rhoS * default_physical.g * Ro^2 / default_physical.sigmaS;

        if default_numerical.version ~= version || ...
           default_numerical.harmonics_qtt ~= harmonics_qtt || ...
           min((data.We - We).^2 + (data.Oh - Oh).^2 + (data.Bo - Bo).^2) < 1e-10    %ismember([We, Oh, Bo], table2array(data(:, {'We', 'Oh', 'Bo'})), 'rows')
            continue;
        end

        bounce = is_there_bounce(recorded_conditions, length_unit);
        % Adding results to table
        data(ii, :) = {bounce, We, Oh, Bo, files_folder(ii).name};
    end
    data = data(~isnan(data.bounce), :);
    writetable(data, sprintf('../2_output/bouncedataharm_qtt=%dversion%d.csv', harmonics_qtt, version));

    
end


function [We, flag] = next_guess(guesses, Oh, Bo)
    % Returns the next best guess given the values that we have available
    
    if size(guesses, 2) == 3
        idx = find(and(guesses.Oh == Oh, guesses.Bo == Bo));
        if ~isempty(idx)
            flag = false;
            We = guesses{idx(1), 'We'};
        elseif inpolygon(Oh, Bo, guesses.Oh, guesses.Bo)
            Mdl = fitrgp(guesses, 'We',  'Standardize', 1);
            flag = true;
            We = predict(Mdl, [Oh, Bo]);
        end

    elseif size(guesses, 2) >3
        guesses_new = guesses(and(guesses.Bo == Bo, guesses.Oh == Oh), :);
        bounces = guesses_new.We(guesses_new.bounce == true);
        no_bounces = guesses_new.We(guesses_new.bounce == false);
        minyes = min(bounces);
        maxno  = max(no_bounces); if isempty(maxno); maxno = 0; end
        if isempty(bounces) || (minyes - maxno)/minyes > 1
            flag = true;
            % Train a classifier and predict where the critical weber might be

            %svc = fitcsvm(guesses{:, {'We', 'Oh', 'Bo'}}, guesses.bounce, ...
            %    'KernelFunction', 'RBF', 'KernelScale', 'auto', 'Standardize', true);
            QDA = fitcdiscr(guesses{:, {'We', 'Oh', 'Bo'}}, guesses.bounce, 'DiscrimType', 'quadratic');
            predictor = @(We) predict(QDA, reshape(We, [length(W), 1]));
            We = logspace(-5, 0, 100);
            bounces = predictor(We);
            We = min(We(bounces == 1));
            if isempty(We)
                We = 0.1;
            end

        else
            We = (minyes + maxno)/2;
            if maxno <= minyes
                if (minyes - maxno)/minyes <= 0.05
                    flag = false;
                else
                    flag = true;
                end
            else
                warning('Inconsistent (nonmonotonic) bouncing behaviour. Stopping procedure');
                flag = false;
            end
            
        end
    else
        warning('Couldnt find a next guess')
        flag = true;
        We = (1 + rand()/2)/500;     
    end
end

function bounce = is_there_bounce(recorded_conditions, length_unit)
    % Calculating the drop's south pole
    harmonics_qtt = length(recorded_conditions{1}.deformation_amplitudes);
    drop_south_pole = @(jj) (cos(pi) .* (length_unit + ...
        sum(collectPl(harmonics_qtt, cos(pi)) .* recorded_conditions{jj}.deformation_amplitudes')) + ...
        recorded_conditions{jj}.center_of_mass) / length_unit >= 0.02;

    % Checking if bounce is present
    if any(arrayfun(drop_south_pole, length(recorded_conditions):-1:1))
        bounce = 1; % True = bounce was present
    else
        bounce = 0; % False = bounce was not present
    end
end