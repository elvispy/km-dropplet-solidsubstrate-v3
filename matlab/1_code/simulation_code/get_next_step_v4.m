function [probable_next_conditions, errortan] = ...
    get_next_step_v4(previous_conditions, dt, contact_points, PROBLEM_CONSTANTS)

    % This function tries to minimize the objetive function by Newton
    % Method
    
    nb_harmonics = previous_conditions{end}.nb_harmonics;
    settings.Fr = PROBLEM_CONSTANTS.froude_nb;
    settings.weights = ones(1, 2*PROBLEM_CONSTANTS.angles_qtt+nb_harmonics)';%[1./(1:nb_harmonics-1) 1./(1:nb_harmonics-1) ...
    settings.M = PROBLEM_CONSTANTS.angles_qtt;
    settings.theta_vector = PROBLEM_CONSTANTS.theta_vector; theta_vector = settings.theta_vector;
    settings.legendre_matrix = PROBLEM_CONSTANTS.precomputed_integrals;
    settings.Oh = PROBLEM_CONSTANTS.Oh;
    %1./(1:nb_harmonics) 1 1 1]';
    
    
    ftm = @(Xn) function_to_minimize(Xn, previous_conditions, dt, contact_points, settings);
    jaccalc = @(Xn) JacobianCalculator(Xn, previous_conditions, dt, contact_points, settings);

    % Initial conditions
    
    Xn = [previous_conditions{end}.deformation_amplitudes(2:end)'; previous_conditions{end}.deformation_velocities(2:end)'; ...
                previous_conditions{end}.pressure_amplitudes(:); previous_conditions{end}.center_of_mass; ...
                previous_conditions{end}.center_of_mass_velocity]; % Initial guess for Newton-Raphson method
    initial_condition = Xn;
    if contact_points >= 0
        % Newton Method
        for m = 1:100
            %perturb = jaccalc(Xn)\ftm(Xn);
            [perturb, ~] =  lsqr(jaccalc(Xn), ftm(Xn), [], 500, [], [], Xn);
            Xnp1 = Xn - perturb;
            %if PROBLEM_CONSTANTS.DEBUG_FLAG == true; plot_condition(2, [0; Xnp1(1:(nb_harmonics-1))]); end
            Xn = Xnp1;
            if norm(ftm(Xnp1)) < 1e-7 || norm(perturb) < 1e-7
                break
            end
        end % end for
        if m == 101
            disp("Newton didn't converge");
            pause
        end
    end % end outer if
    
    if norm(initial_condition-Xn)/norm(initial_condition) >= 1e-2/nb_harmonics
        warning('Finished newton method. The solution is too far away from the previous condition'); 
    end
    
    % Lets try with the same pressure
    probable_next_conditions = previous_conditions{end};
    M = probable_next_conditions.nb_harmonics;
    probable_next_conditions.deformation_amplitudes = [0; Xn(1:(M-1))]';
    probable_next_conditions.deformation_velocities = [0; Xn(M:(2*M-2))]';
    probable_next_conditions.pressure_amplitudes    = Xn((end-2-M):(end-2))';
    probable_next_conditions.center_of_mass         = Xn(end-1);
    probable_next_conditions.center_of_mass_velocity= Xn(end);
    probable_next_conditions.dt = dt;
    probable_next_conditions.current_time = previous_conditions{end}.current_time + dt;
    probable_next_conditions.contact_points = contact_points;
%     probable_next_conditions.contact_radius = ...
%         round(sin(PROBLEM_CONSTANTS.theta_vector(contact_points)) * (1 + ...
%         sum(arrayfun(@(idx) (-1)^(idx) * Xnp1(idx), 1:length(Xnp1)))), 10);
    z_calculator = @(angle) cos(angle) .* (1 + sum(probable_next_conditions.deformation_amplitudes' .* ...
        collectPl(nb_harmonics, cos(angle)))) + probable_next_conditions.center_of_mass;
    errortan = 0;
    angles_check = theta_vector(theta_vector > pi/2);
    if contact_points > 0
        errortan = (z_calculator(theta_vector(contact_points+1)) - ...
            z_calculator(theta_vector(contact_points)))/...
            (theta_vector(contact_points) - theta_vector(contact_points+1));
        angles_check = angles_check(angles_check < theta_vector(contact_points));
    end
   
    check = any(z_calculator(angles_check) < 0);
    
    if (check || contact_points < 0); errortan = inf; end
    
    
end % end main function definition




