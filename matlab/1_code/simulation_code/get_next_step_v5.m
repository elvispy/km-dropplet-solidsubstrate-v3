function [probable_next_conditions, errortan, mat_inverses] = ...
    get_next_step_v5(previous_conditions, dt, contact_points, PROBLEM_CONSTANTS)

    % This function tries to minimize the objetive function by Newton
    % Method
    
    nb_harmonics = previous_conditions{end}.nb_harmonics;
    settings.Fr = PROBLEM_CONSTANTS.froude_nb;
    settings.weights = ones(1, 2*PROBLEM_CONSTANTS.angles_qtt+nb_harmonics)';%[1./(1:nb_harmonics-1) 1./(1:nb_harmonics-1) ...
    settings.M = PROBLEM_CONSTANTS.angles_qtt;
    settings.theta_vector = PROBLEM_CONSTANTS.theta_vector; theta_vector = settings.theta_vector;
    settings.legendre_matrix = PROBLEM_CONSTANTS.precomputed_integrals;
    settings.Oh = PROBLEM_CONSTANTS.Oh;
    mat_inverses.version = PROBLEM_CONSTANTS.version;
    %1./(1:nb_harmonics) 1 1 1]';
    
    
    ftm = @(Xn) PROBLEM_CONSTANTS.function_to_minimize(Xn, previous_conditions, dt, contact_points, settings);
    jaccalc = @(Xn) PROBLEM_CONSTANTS.jacobian_calculator(Xn, previous_conditions, dt, contact_points, settings);

    % Initial conditions
    
    Xn = [previous_conditions{end}.deformation_amplitudes(2:end)'; previous_conditions{end}.deformation_velocities(2:end)'; ...
                previous_conditions{end}.pressure_amplitudes(:); previous_conditions{end}.center_of_mass; ...
                previous_conditions{end}.center_of_mass_velocity]; % Initial guess for Newton-Raphson method
    initial_condition = Xn;
    best_solution = Xn;
    
    %Start newton method only if contact_points >= 0 (discarding degenerate cases)
    if contact_points >= 0
        best_value = norm(ftm(Xn));
        % Newton Method
        for m = 1:100
            %perturb = jaccalc(Xn)\ftm(Xn);
            [perturb, mat_inverses] = my_lsqr(jaccalc, ftm, Xn, mat_inverses, dt);
            Xnp1 = Xn - perturb;
            %if PROBLEM_CONSTANTS.DEBUG_FLAG == true; plot_condition(2, [0; Xnp1(1:(nb_harmonics-1))]); end
            Xn = Xnp1;
            score = norm(ftm(Xnp1));
            if score < best_value
                best_value = score;
                best_solution = Xnp1;
            end
            if (score < 1e-10 || norm(perturb) < 1e-10)
                break
            end
        end % end for
    end % end outer if
    
    if norm(initial_condition-best_solution)/norm(initial_condition) >= 1/nb_harmonics
        %warning('Finished newton method. The solution is too far away from the previous condition'); 
    end
    
    % Lets try with the same pressure
    probable_next_conditions = previous_conditions{end};
    M = probable_next_conditions.nb_harmonics;
    probable_next_conditions.deformation_amplitudes = [0; best_solution(1:(M-1))]';
    probable_next_conditions.deformation_velocities = [0; best_solution(M:(2*M-2))]';
    probable_next_conditions.pressure_amplitudes    = best_solution((end-2-M):(end-2))';
    probable_next_conditions.center_of_mass         = best_solution(end-1);
    probable_next_conditions.center_of_mass_velocity= best_solution(end);
    probable_next_conditions.dt = dt;
    probable_next_conditions.current_time = previous_conditions{end}.current_time + dt;
    probable_next_conditions.contact_points = contact_points;

    z_calculator = @(angle) -cos(angle) .* (1 + sum(probable_next_conditions.deformation_amplitudes' .* ...
        collectPl(nb_harmonics, cos(angle)))) + probable_next_conditions.center_of_mass; % minus sign = south pole reference
    errortan = 0;

    angles_check = theta_vector(theta_vector < pi/2);
    if contact_points > 0
        errortan = (z_calculator(theta_vector(contact_points+1)) - ...
            z_calculator(theta_vector(contact_points)));%/...
            %(theta_vector(contact_points+1) - theta_vector(contact_points));
        angles_check = angles_check(angles_check > theta_vector(contact_points));
    end
   
    check = any(z_calculator(angles_check) < 0);
    
    if (check || contact_points < 0); errortan = inf; end
    
    
end % end main function definition

function [perturb, mat_inverse] = my_lsqr(jaccalc, ftm, Xn, mat_inverse, dt)
    A = jaccalc(Xn);
    fieldName = strrep(strrep(sprintf('dt%.2e', dt), '-', ""), ".", "");
    S = size(A);
    %if S(1) == S(2) % If it's square matrix, just invert it
        % If there's already one in PROBLEM_CONSTANTS, use it
    if isfield(mat_inverse, fieldName)
        perturb = mat_inverse.(fieldName) * ftm(Xn);
    else
        if mat_inverse.version == 3 && S(1) == S(2) % Only take inverse if it's a linearized model
            mat_inverse.(fieldName) = inv(A);
        end
        perturb = A\ftm(Xn);
    end
    % else
    %     [perturb, ~] = lsqr(A, ftm(Xn), [], 500, [], [], Xn);
    % end
end