% Author: Elvis Aguero
% email: elvisavfc65@gmail.com
% Date: January 13th, 2023.

function solve_motion_v2(varargin)
    
    %
    %Tries to solve the full kinematic match between a dropplet
    % and a solid substrate in vacuum conditions.
    % Arguments:
    % varargin must contain at most three elements, and all of them must be structs whose fields may be left
    % empty or assume nan values to adopt default values, when possible.
    %  - varargin{1} = Physical parameters for the simulation. All dimensions are in cgs units. Can contain the following fields:
    %                                (dim, default) <-- Default shown as ? if variable is mandatory
    %      1) undisturbed_radius     (cm, 1)          = radius of the undisturbed sphere 
    %      2) com_height             (cm, contact)    = Initial height of the COM of the drop. Default is iminent contact
    %      3) impact_velocity        (cm/s, ?)        = Initial velocity of the COM (negative = towards impact)
    %      4) amplitudes             (cm, all zeroes) = Spectral legendre amplitudes for the surface of the drop. 1-indexed
    %      5) amplitudes_vel         (cm, all zeroes) = Spectral legendre velocities for the surface of the drop. 1-indexed
    %      6) pressure_amplitudes    (. , all zeroes) = Initial pressure amplitudes in spectral coordinates. 0-indexed
    %      7) initial_contact_points (adim, 0)        = Initial number of contact points
    %      8) rhoS                   (kg/cm^3, 0.988) = The density of the fluid inside the droplet. Default = water
    %      9) sigmaS                 (., 72.29)       = Surface tension of the fluid inside the droplet. Default = water
    %      10) g                     (cm/s^2, 9.81e+2)= Gravitational constant.
    %  - varargin{2} = Numerical parameters for the simulation
    %      1) harmonics_qtt          (adim, ?)        = Number of spectral amplitudes that describe the motion. 
    %      2) angular_sampling (adim, harmonics_qtt+1)= Number of angles that describe the shape of the drop.
    %      3) simul_time             (s, 15e-3)       = Maximum simulation time allowed.
    %  - varargin{3} = Other options. 
    %      1) live_plotting          (bool, false)    = whether or not to plot real-time results (more consuming)
    %      2) debug_flag             (bool, false)    = Verbose real-time info for the simulation (experimental feature)

    if nargin >= 3
        default_options = struct('live_plotting', false, 'debug_flag', false);
        A = fieldnames(varargin{3})
        for ii = 1:length(A)
            default_options.(A{ii}) = varargin{3}.(A{ii})
        end
    end
    if nargin >= 2
        default_numerical = struct('simul_time', 15e-3, 'harmonics_qtt', nan, 'angular_sampling', nan);
        if isstruct(varargin{2}) == false; error('Numerical values is not a struct'); end
        A = fieldnames(varargin{2})
        for ii = 1:length(A)
            default_numerical.(A{ii}) = varargin{2}.(A{ii})
        end
        if length(A) < 3; default_numerical.angular_sampling = default_numerical.harmonics_qtt + 1; end 
    end
    if nargin >= 1

    %% Handling default arguments. All units are in cgs.
    
    undisturbed_radius = .1;  % Radius of the undeformed spherical sphere 
    initial_height = Inf;    % Initial position of the sphere center of mass of the sphere (Inf = start at imminent contact)
    initial_velocity = -10; % Initial velocity of the sphere in cm/s
    initial_amplitudes = Inf; % Initial amplitudes of the dropplet (Default = undisturbed) OBS: First index is A_1, not A0
    initial_contact_points = 0;
    amplitudes_velocities = [];
    rhoS = 0.998;            % Sphere's density
    sigmaS = 72.20;          % Sphere's Surface Tension
    g = 9.8065e+2;           % Gravitational constant
    harmonics_qtt = 15;      % Number of harmonics to be used 
    angular_sampling = harmonics_qtt + 1; 
    max_dt = 0;         % maximum allowed temporal time step
    %angle_tol =  pi * 5/harmonics_qtt;
    simulation_time = 15e-3; % Maximum allowed total time in seconds
    live_plotting = true; % Whether to plot or not the live results

    % Dimensionless Units
    length_unit = undisturbed_radius;
    time_unit = sqrt(rhoS * length_unit^3 / sigmaS); %undisturbed_radius/velocity_unit; % Temporal dimensionless number
    velocity_unit = length_unit/time_unit; % abs(initial_velocity);    
    pressure_unit = rhoS * velocity_unit^2;
    froude_nb   = length_unit/(g*time_unit.^2);
    weber_nb    = rhoS * undisturbed_radius.^3/ (sigmaS*time_unit.^2); % Weber's number of the dropplet
    mS = rhoS * length_unit^3;
    mass_unit = mS;
    
    f = @(n)  sqrt(n .* (n+2) .* (n-1));
    omegas_frequencies = f(1:harmonics_qtt)';
    f2 = @(n) (1 - n ./ ((angular_sampling+n):-1:(1+n))) * pi * (angular_sampling + n)/angular_sampling;
    syms x;
    
    % we choose the angles to be the zeros of the last legendre Polynomial
    theta_vector = [pi; acos(double(vpasolve(legendreP(angular_sampling-1, x))))]'; clear x;
    %theta_vector = f2(harmonics_qtt); % The bigger the argument of f2, the more uniform is the distribution.
    % The smallest the argument, the more skewed is towards pi. For an uniform distribution, use linspace(pi, 0, angular_sampling);

    
    % % Initial conditions
    % Set dropplet's sphere height initial conditions
    get_initial_height = @(amplitudes) 1 - sum(arrayfun(@(idx) amplitudes(idx) * (-1.0)^(idx), 1:length(amplitudes)));
    
    if length(initial_amplitudes) ~= harmonics_qtt
        initial_amplitudes = zeros(1, harmonics_qtt);
    else
        initial_amplitudes = initial_amplitudes/length_unit;
    end
    if initial_height == Inf
        initial_height = get_initial_height(initial_amplitudes);
    else
        assert(initial_height >= get_initial_height(initial_amplitudes));
        initial_height = initial_height/length_unit;
    end
    
    initial_velocity_adim = initial_velocity/velocity_unit;

    if length(amplitudes_velocities) ~= harmonics_qtt
        amplitudes_velocities = zeros(1, harmonics_qtt);
    else
        amplitudes_velocities = amplitudes_velocities/velocity_unit;
    end

    % tan_tol = tan(min_angle);

    initial_pressure_coefficients = zeros(1, harmonics_qtt+1) / pressure_unit; % Just to emphasize the units of these coefficients.
    contact_points = 0;
    if max_dt == 0
        dt = round(time_unit/(20 * harmonics_qtt^(1/2)), 1, 'significant')/time_unit; % /(10 * harmonics_qtt^(1/2))
    else 
        dt = max_dt/time_unit; 
    end
    
    
    initial_time = 0;
    current_time = initial_time/time_unit;
    final_time = simulation_time/time_unit;
    current_index = 2;%  This integer points to the next available index in variables that are going to 
                      %  export data (index 1 is for initial conditions)
    maximum_index = ceil((final_time - initial_time)/dt) + 4;
    number_of_extra_indexes = 0;

    %contact_points = 0;%  Initial number of contact points
    contact_time = 0;% TODO: Lab contact time and so on?
    labcTime = 0;%To record contact time
    maxDef = 0; %to record max deflection;
    firstFallFlag = true; %Boolean to record maximum deflection time
    contactFlag = false; %To record first contact time
    velocityOutRecorded = false; labvelocityOutRecorded = false;  % To check if Em_out was recorded
    grow_dt = false;%  THis variable controls how fast dt can grow
    iii = 0; jjj = 0;%  Indexes to keep track how small is dt compared to max_dt

      
    legendre_matrix = precompute_integrals(theta_vector, harmonics_qtt);
    % Constants of the problem formulation
    PROBLEM_CONSTANTS = struct("froude_nb", froude_nb, "weber_nb", weber_nb, ...
        "nb_harmonics", harmonics_qtt, ...
        "omegas_frequencies", omegas_frequencies, ...
        "angles_qtt", harmonics_qtt + 1, ... % number of angles 
        "pressure_unit", pressure_unit, ...
        "theta_vector", theta_vector, ... % set of fixed angle vectors
        "precomputed_integrals", legendre_matrix, ... % integral os legendre polynomials in angle intervals
        "function_to_minimize", @function_to_minimize_v3, ... % v1 = fully nonlinear integration on disk, v2 = nonlinear with spherical approximation, v3 = linearised version of v2
        "jacobian_calculator", @JacobianCalculator_v3, ... % has to be the same version as function to minimize
        "DEBUG_FLAG", true); %true = plot and save video

    %current_conditions = cell(1, 1); % probable_next_conditions = Vector{ProblemConditions}(undef, 5);
    current_conditions = ProblemConditions_v2( ...
        harmonics_qtt, ...
        initial_amplitudes, ...
        amplitudes_velocities, ...
        initial_pressure_coefficients, ...
        current_time, ...
        dt, ...
        initial_height, ...
        initial_velocity_adim, initial_contact_points); % Last argument is contact radius
 
    previous_conditions = {current_conditions, current_conditions}; 
    % TODO: Define this array properly to implement BDF2.
    previous_conditions{1}.current_time = previous_conditions{2}.current_time - dt;
    previous_conditions{1}.center_of_mass_velocity = ...
        previous_conditions{2}.center_of_mass_velocity + dt/froude_nb;
    previous_conditions{1}.center_of_mass = ...
        previous_conditions{2}.center_of_mass - previous_conditions{2}.center_of_mass_velocity * dt;
    
    g = @(t, idx) current_conditions.deformation_amplitudes(idx) * cos(f(idx) * t) ...
        + current_conditions.deformation_velocities(idx)/(f(idx)+1e-30) * sin(f(idx) * t); 

    for idx = 1:harmonics_qtt
        previous_conditions{1}.deformation_amplitudes(idx) = g(-dt, idx);
        previous_conditions{1}.deformation_velocities(idx) = (g(0, idx) - g(-2*dt/1000, idx))/(2*dt/1000);
    end
    % 1-st order set
    previous_conditions = previous_conditions(end);

    % % Preparing post-processing
    % TODO: Write post processing variables

   %  Preallocate variables that will be exported (All of them have units!)
   recorded_conditions =cell(maximum_index, 1); % Vector{ProblemConditions}(undef, (maximum_index, )); 
   give_dimensions_v2 = @(X) ProblemConditions_v2( ...
       X.nb_harmonics, ...
        X.deformation_amplitudes * length_unit, ...
        X.deformation_velocities * velocity_unit, ...
        X.pressure_amplitudes * (mass_unit * length_unit / (time_unit^2 * length_unit^2)), ...
        X.current_time * time_unit, ...
        X.dt * time_unit, ...
        X.center_of_mass * length_unit, ...
        X.center_of_mass_velocity * velocity_unit, ...
        X.contact_points); 

     recorded_conditions{1} = give_dimensions_v2(previous_conditions{end});
     recorded_times = zeros(maximum_index, 1); recorded_times(1) =  current_time * time_unit;
%     
%    %  Coefficient of restitution
%     mechanical_energy_in = NaN;
%     mechanical_energy_out = NaN; % TODO: Lab COef of restitution?

    indexes_to_save = zeros(maximum_index, 1); indexes_to_save(1) = 1;
    current_to_save = 2;
    
    if PROBLEM_CONSTANTS.DEBUG_FLAG == true
        file_name2 = 'DropAndSubstrate.mp4';
        vidObj = VideoWriter(file_name2,'MPEG-4');
        set(vidObj,'FrameRate',10)
        open(vidObj);
    end
    
    

    
    % Define advancing step function: 
    %   v4 = classic newton method
    %   v5 = Newton method with best value estimator
    advance_one_step = @(a, b, c, d) get_next_step_v5(a, b, c, d);
    %% Starting main loop
    while ( current_time < final_time) 
        % First, we try to solve with the same number of contact points
        probableNextConditions = cell(5, 1);
        errortan = Inf * ones(1, 5);
        recalculate = false;
        
        [probableNextConditions{3}, errortan(3)] = advance_one_step(previous_conditions, dt, ...
            contact_points, PROBLEM_CONSTANTS);
        [probableNextConditions{4}, errortan(4)] = advance_one_step(previous_conditions, dt, ...
            contact_points+1, PROBLEM_CONSTANTS);
        [probableNextConditions{2}, errortan(2)] = advance_one_step(previous_conditions, dt, ...
            contact_points-1, PROBLEM_CONSTANTS);
            
        if (abs(errortan(3)) > abs(errortan(4)) || abs(errortan(3)) > abs(errortan(2)))
            if abs(errortan(4)) <= abs(errortan(2))
                %Now lets check with one more point to be sure
                [~, errortan(5)] = advance_one_step(previous_conditions, dt, ...
        contact_points+2, PROBLEM_CONSTANTS);

                if abs(errortan(4)) < abs(errortan(5))
                    %Accept new data
                    
                    current_conditions = probableNextConditions{4};
                    contact_points = contact_points + 1;
                else
                    %time step too big
                    recalculate = true;
                end
            else
                %now lets check if errortan is good enough with one point
                %less
                [~, errortan(1)] = advance_one_step(previous_conditions, dt, ...
        contact_points-2, PROBLEM_CONSTANTS);


                if abs(errortan(2)) < abs(errortan(1))
                    %Accept new data
                    current_conditions = probableNextConditions{2};
                    contact_points = contact_points - 1;
                else
                    recalculate = true;
                end
            end %End of (errortan(4) < errortan(2))

        else %the same number of contact points is best    
            if errortan(3) == Inf % ALl errors are infinity
                recalculate = true;
            else
                %Accept new data
                current_conditions = probableNextConditions{3};
                
            end
        end %
        
        if recalculate == true
            dt = dt/2;
            disp("Se dividio dt por 2");
            % Refine time step in index notation 
            iii = iii + 1; jjj = 2 * jjj;
        else
            previous_conditions = {previous_conditions{2:end} current_conditions};
            
            current_time = current_time + dt; jjj = jjj + 1;
            if mod(jjj, 2) == 0 && grow_dt == true
                jjj = floor(jjj/2); 
                iii = iii - 1;
                % Increase time step
                dt = 2 * dt;
                % Decrease the number of time you can make dt bigger
                grow_dt = false;
            end

            %  TODO: Update Indexes if necessary

            % TODO: % Stored data
            recorded_conditions{current_index} = give_dimensions_v2(current_conditions);
            recorded_times(current_index) = current_time * time_unit;
            current_index = current_index + 1; % Point to the next space in memory 

            % If we are in a multiple of max_dt, reset indexes
            if jjj == 2^iii
                jjj = 0;
                grow_dt = true;
                indexes_to_save(current_to_save) = current_index - 1;
                current_to_save = current_to_save + 1;
            else
                number_of_extra_indexes = number_of_extra_indexes + 1;
            end

            if live_plotting == true
                % Do some live plotting here

                plot_title = sprintf(" t = %-8.5f (ms), Contact points = %-2g deg, \n v_k = %-8.5f cm/s, z_k = %-8.5f cm\n", ...
                   1e+3 * current_time * time_unit, current_conditions.contact_points, ...
                        current_conditions.center_of_mass_velocity * velocity_unit, ...
                        current_conditions.center_of_mass* length_unit);
                h = plot_condition(1, current_conditions, 1.25, plot_title, PROBLEM_CONSTANTS.theta_vector);
                
                if PROBLEM_CONSTANTS.DEBUG_FLAG == true
                    currFrame = getframe(h);
                    writeVideo(vidObj,currFrame);
                end

            else
                % Do some real-time variable updating here
            end
            
            %%%%%%%%%%%%
            %%%%%%%%%%%%
            % ANALISE CONTACT TIME, MAXIMUM DEFLECTION, COEF OF RESTITUTION
            if (contactFlag == false && current_conditions.contact_points > 0) % if contact began, start counting
                contactFlag = true;
                
                impact_velocity = round(recorded_v(current_index - 2, 1), 10); % Store initial impact velocity
                zeroPotential = lenght_unit; % Store our zero-Potential 
                Em_in = 1/2 * mS * ((recorded_conditions{current_index - 2}.center_of_mass_velocity)^2); % Mechanical Energy In;

            elseif (contactFlag == true && ...
                    (velocityOutRecorded == false || labvelocityOutRecorded == false)) % record last contact time
                if recorded_conditions{current_index - 1}.center_of_mass > lenght_unit 
                    if labvelocityOutRecorded == false
                        labEm_out = 1/2 * mS * ((recorded_conditions{current_index - 2}.center_of_mass_velocity)^2) ...
                            + mS*g*(recorded_conditions{current_index - 2}.center_of_mass - zeroPotential); % Mechanical Energy Out;
                        labvelocityOutRecorded = true;
                    end
                elseif labvelocityOutRecorded == false
                    labcTime = labcTime + dt;
                end
                if current_conditions.contact_points == 0 
                    if velocityOutRecorded == false
                        Em_out = 1/2 * mS * ((recorded_conditions{current_index - 2}.center_of_mass_velocity)^2) ...
                            + mS*g*(recorded_conditions{current_index - 2}.center_of_mass - zeroPotential); % Mechanical Energy Out;
                        velocityOutRecorded = true;
                    end
                elseif velocityOutRecorded == false
                    contact_time = contact_time + dt;
                end
            end

            if (contactFlag == true && recorded_conditions{current_index - 1}.center_of_mass_velocity >= 0) % Record maximum deflection
                if (firstFallFlag == true) %If velocity has changed sign, record maximum 
                    %deflection only once
                    maxDef = maxDef - recorded_conditions{current_index - 2}.center_of_mass;
                    firstFallFlag = false;
                end
            end
            if  recorded_conditions{current_index - 1}.center_of_mass_velocity < 0 && firstFallFlag == false 
                if velocityOutRecorded    == false; contact_time    = inf; end
                if labvelocityOutRecorded == false; labcTime = inf; end
                
            end

        end

    end
    
    if PROBLEM_CONSTANTS.DEBUG_FLAG == true
        close(vidObj);
    end
    
    
    file_path = fullfile("../0_data/", 'manual');
    if exist(file_path, 'dir') ~= 7 % CHeck if folder simulations exists
        mkdir(file_path); % If not, create it
    end
    file_name = fullfile(file_path, sprintf("simulation %s.mat", datestr(now, 0)));
    % Post processing
    impact_parameters = struct('
    indexes_to_save = indexes_to_save(1:(current_to_save-1));
    recorded_conditions = recorded_conditions(indexes_to_save);
    save(file_name, 'recorded_conditions', 'recorded_times', 'length_unit', ...
        'velocity_unit', 'pressure_unit', 'time_unit', 'froude_nb', 'weber_nb', 'PROBLEM_CONSTANTS');

 
    
end


        
