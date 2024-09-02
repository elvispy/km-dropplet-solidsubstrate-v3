function res = JacobianCalculator_v1(Xn, previous_conditions, dt, contact_points, settings)

    %Ps = settings.Ps;
    %Flat = settings.Flat;
    Oh = settings.Oh;
    weights = settings.weights; % weights of all the equations present in the system
    n = length(previous_conditions); % Determines the order of the method
    theta_vector = settings.theta_vector; if size(theta_vector, 2) > 1; theta_vector = theta_vector'; end
    N = length(theta_vector);% Number of contact angles in sphere
    if n > 2 || n < 1; throw("Hey!"); end

    M = previous_conditions{end}.nb_harmonics;
    if n == 1
        coefs = [-1.0, 1.0];
    elseif n == 2
        rk = dt/previous_conditions{end}.dt;
        ak = (1+2*rk)/(1+rk);
        bk = -(1+rk);
        ck = rk^2/(1+rk);
        coefs = [ck, bk, ak]; 
    end

    %X = [Al, \dot{Al}, Bl, z, \dot{z}]

    % We proceed by newton-raphson

    % FIRST ROW BLOCK (Deformation evolution pt1)
    R1 = [coefs(end) * eye(M-1), -dt * eye(M-1), zeros(M-1, M+1), zeros(M-1, 2)];

    % Second ROW BLOCK (Deformation evolution pt2 )
    Aidx = 2:M;
    %Bidx = 0:N;
    D1N = diag(Aidx .* (Aidx  + 2) .* (Aidx - 1));
    D2N = [zeros(M-1, 2), diag(Aidx)]; % B0 and B1 do not contribute to the motion of the drop
    R2 = [dt*D1N, diag(coefs(end)+2*dt*Oh*(2*Aidx+1) .* (Aidx-1)), dt * D2N, zeros(M-1, 2)];


        
        % third ROW BLOCK (flat surface on contact angle condition)
        %L = 10;
        %theta_i2 = linspace(theta_contact, pi, Flat)'; x_i2 = cos(theta_i2);
        %P2 = sin(theta_i2) .* (collectPl(N, x_i2)' + x_i2 .* collectdnPl(N, x_i2)');
        theta_contact = theta_vector(1:contact_points);
        P2 = collectPl(M, cos(theta_contact))';
        B2 = cos(theta_contact) .* P2(:, 2:end); % Discard A1
        R3 = [B2, zeros(contact_points, M-1), zeros(contact_points, M+1), ones(contact_points, 1), zeros(contact_points, 1)];
        %R3 = R3(1:contact_points, :);
    
        % Fourth ROW BLOCK (No pressure outside condition)
        theta_i = theta_vector((contact_points+1):end);
        cosines = cos(theta_i);
        P1 = collectPl(M, cosines)';
        B1 = [ones(size(cosines)), P1];
        R4 = [zeros(N-contact_points, 2*(M-1)), B1, zeros(N-contact_points, 2)]/1e+3;


        % Fifth ROW BLOCK (Center of mass condition)
        %oneN = (-1).^(Aidx);
        %R5 = [-oneN, zeros(1, 2*N), 1, 0];

        % Sixth ROW BLOCK (Center of mass diff equation)
        R6 = [zeros(1, M-1), zeros(1, M-1), zeros(1, M+1), coefs(end), -dt];
        
        % Seventh ROW BLOCK (\dot{v} = integral of pressure)
        integral_coefficients = sum(settings.legendre_matrix(1:contact_points, :), 1);
        current_pressures = Xn((end-M-2):(end-2))';
        current_heights     = Xn(end-1);
        Zl = 3/2 * current_heights.^2 * integral_coefficients;
        w = 3 * current_heights * sum(integral_coefficients .* current_pressures);
        R7 = [zeros(1, M-1), zeros(1, M-1), dt * Zl, dt*w, coefs(end)];
        
        % Seventh ROW BLOCK (center of mass velocity diff equation)
        %oneN = (-1).^(Aidx);
        
        %A = (1 + sum(current_deformation .* oneN))^2;
        %Cl = manual_intPnxm3(N, -1, cos(theta_contact))';
        %B = sum(current_pressures .* Cl);

        %R7 = [-dt * 3 * (current_deformation + oneN) * B, ...
        %    zeros(1, N-1), -dt*(3/2 * A) * Cl, 0, coefs(end)];
%     else
%         R3 = [zeros(N+1, 2*(N-1)), eye(N+1), zeros(N+1, 2)];
%         R4 = zeros(0, 3*N+1);
%         R5 = R4;
%         R6 = [zeros(1, 3*N-1), coefs(end), -dt];
%         R7 = [zeros(1, 3*N-1),   0, coefs(end)];
%     end
        
    res = [R1;R2;R3;R4;R6;R7];
    
end