function res = function_to_minimize_v2(Xn, previous_conditions, dt, contact_points, settings)
    % Returns the function to be minimized for the Newton-Raphson Method
    
    Fr = settings.Fr;
    %Ps = settings.Ps;
    %Flat = settings.Flat;
    %M = settings.M; % Number of contact angles in sphere
    %weights = settings.weights;
    n = length(previous_conditions); % Determines the order of the method
    theta_vector = settings.theta_vector; if size(theta_vector, 2) > 1; theta_vector = theta_vector'; end
    %N = length(theta_vector);
    if n > 2 || n < 1; throw("Hey!"); end

    M = previous_conditions{end}.nb_harmonics;
    %dt = Xn.dt;
    if n == 1
        coefs = [-1.0, 1.0];
    elseif n == 2
        rk = dt/previous_conditions{end}.dt;
        ak = (1+2*rk)/(1+rk);
        bk = -(1+rk);
        ck = rk^2/(1+rk);
        coefs = [ck, bk, ak]; 
    end
    
    extract_symbol = @(jj, field) previous_conditions{jj}.(field);
    previous_velocities = zeros(M, n);
    previous_deformation= zeros(M, n);
    previous_pressures  = zeros(M+1, n);
    previous_COM        = zeros(1, n);
    previous_COM_vel    = zeros(1, n);
    for ii = 1:n
        previous_velocities(:, ii) = reshape(extract_symbol(ii, 'deformation_velocities'), M, 1);
        previous_pressures(:, ii)  = reshape(extract_symbol(ii, 'pressure_amplitudes'), M+1, 1);
        previous_deformation(:,ii) = reshape(extract_symbol(ii, 'deformation_amplitudes'), M, 1);
        previous_COM(1, ii)        = extract_symbol(ii, 'center_of_mass');
        previous_COM_vel(1, ii)    = extract_symbol(ii, 'center_of_mass_velocity');
    end
    
    idxs = 1:(M-1);
    previous_velocities  = previous_velocities(2:end, :);
    previous_deformation = previous_deformation(2:end, :);
    
    current_deformation  = Xn(idxs);
    current_velocities   = Xn(idxs + (M-1));
    current_pressures    = Xn((end-2-M):(end-2));
    center_of_mass       = Xn(end-1);
    center_of_mass_vel   = Xn(end);
    
    % FIRST ROW BLOCK (Deformation evolution pt1)
    R1 = (sum(coefs .* [previous_deformation, ...
        current_deformation], 2) - dt * current_velocities);

    % Second ROW BLOCK (Deformation evolution pt2 )
    Aidx = (2:M)';
    
    D1N = Aidx .* (Aidx + 2) .* (Aidx - 1);
    D1N2 = 2*settings.Oh*(Aidx-1) .* (2*Aidx+1);
    
    R2 =  (sum(coefs .* [previous_velocities, current_velocities], 2) ...
         + dt * (Aidx .* current_pressures(3:end) ...
         + D1N2 .* current_velocities ... % Added viscocity term
         + D1N .* current_deformation));
     
    %if theta_contact < pi
        % Third ROW BLOCK (flat surface on contact angle condition)
        %L = 10; % This, too, must be unified.
        theta_contact = theta_vector(1:contact_points); %theta_i2 = reshape(linspace(theta_contact, pi, Flat), Flat, 1); 
        
        cosines = cos(theta_contact);
        P2 = collectPl(M, cosines)';
        P2 = P2(:, 2:end); % Discard A1
        R3 = cosines .* (1 + sum(P2 * current_deformation, 2)) + center_of_mass;

    
        % fourth ROW BLOCK (No pressure outside condition)
        theta_i = theta_vector((contact_points+1):end);
        P1 = collectPl(M, cos(theta_i))';
        P1 = [ones(size(theta_i)), P1];
        R4 = sum(P1 * current_pressures, 2);

        
        % Fifth ROW BLOCK (Center of mass condition)
        %oneN = (-1).^(Aidx);
        %alternating_sum = sum(oneN .* current_deformation);
        %R5 =  (center_of_mass - alternating_sum - 1);

        % Sixth ROW BLOCK (Center of mass diff equation)
        R6 =  sum(coefs .* [previous_COM, center_of_mass] , 2) - dt * center_of_mass_vel; %[zeros(1, 3*N-1), coefs(end), -dt];

        % Seventh ROW BLOCK (center of mass velocity diff equation, exact pressure integral on contact area)
        
        current_pressures = Xn((end-M-2):(end-2))';
        idx = 2:M;
        Cl = 3 * idx .* (idx-1) ./ (2.*idx - 1) ./ (2.*idx + 1);
        Dl = 3 * (idx + 2) .* (idx + 1) ./ (2 * idx + 3) ./ (2.*idx + 1);
        Xl = (Cl .* current_pressures(2:M) - Dl .* [current_pressures(4:(M+1)), 0])';
        %Yl = [0, 1, [Cl(2:end) .* current_amplitudes(3:M), 0] - [0, Dl(1:(M-2)) .* current_amplitudes(2:(M-1))]];
        
        
        R7 =  (sum(coefs .* [previous_COM_vel, center_of_mass_vel], 2) - ...
            dt * (-1/Fr - current_pressures(2) + sum(current_deformation .* Xl)));
        
            % [-dt * 3 * (current_deformation + oneN) * B, ...
            %zeros(1, N-1), -dt*(3/2 * A) * Cl, 0, coefs(end)];
%     else
%         R3 = current_pressures;
%         R4 = zeros(0, 1);
%         R5 = zeros(0, 1);
%         R6 = sum(coefs .* [previous_COM, center_of_mass] , 2) - dt * center_of_mass_vel;
%         R7 = sum(coefs .* [previous_COM_vel, center_of_mass_vel], 2) + dt/Fr;
%     end
    %res = weights .* [R1;R2;R3;R4;R6;R7];
    res = [R1;R2;R3;R4;R6;R7];
end