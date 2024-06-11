function tests = function_to_minimize_v1_test
    tests = functiontests(localfunctions);
    %disp("lol")
end

function setupOnce(testCase)
    import matlab.unittest.TestCase
    import matlab.unittest.constraints.IsEqualTo
    import matlab.unittest.constraints.RelativeTolerance
    
    M = 50; dt = 1;
    testCase.TestData.dt = 1;
    testCase.TestData.M = M; % nb_harmonics
    testCase.TestData.N = M+1; % nb_angles
    settings.Fr = 1;
    settings.theta_vector = linspace(pi, 0, testCase.TestData.N);
    settings.legendre_matrix = precompute_integrals(...
        settings.theta_vector, testCase.TestData.M);
    settings.weights = 0;
    testCase.TestData.settings = settings;
    
    testCase.TestData.current_conditions = ProblemConditions_v2( ...
        testCase.TestData.M, ...
        zeros(M, 1), ...
        zeros(M, 1), ...
        zeros(M+1, 1), ...
        0, ...
        dt, ...
        1, ...
        -1, 0); 
    
    previous_conditions = {testCase.TestData.current_conditions, testCase.TestData.current_conditions};
    previous_conditions{1}.current_time = -dt;
    previous_conditions{2}.current_time = - 2*dt;
    
    testCase.TestData.previous_conditions = previous_conditions; 
    
end

function testJacobian(testCase)
    % No contact points, zero initial conditions, zero testing value
    dt = 1;
    contact_points = 0;
    M = testCase.TestData.M;
    
    %current_conditions  = testCase.TestData.current_conditions;
    previous_conditions = testCase.TestData.previous_conditions;
    
    %initial conditions
    Xn = [zeros(M-1, 1); zeros(M-1, 1); zeros(M+1, 1); 1; 1];

    L = length(Xn);
    ftm = @(Xn) function_to_minimize_v1(Xn, previous_conditions, dt, contact_points, testCase.TestData.settings);
    expected_value = zeros(length(ftm(Xn)), length(Xn));
    COL = ftm(Xn);
    ddt = 1e-4;
    for jj = 1:L
        perturb = zeros(L, 1); perturb(jj) = ddt;
        expected_value(:, jj) = (ftm(Xn + perturb) - COL)/ddt;
    end
    actual_value = JacobianCalculator_v1(Xn, previous_conditions, dt, contact_points, testCase.TestData.settings);
    
    verifyEqual(testCase, actual_value, expected_value, 'RelTol', 1e-5); 
end


function testJacobianRandomInitialCondition(testCase)
    % No contact points, zero initial conditions, zero testing value
    dt = 1;
    contact_points = 0;
    M = testCase.TestData.M;
    
    %current_conditions  = testCase.TestData.current_conditions;
    previous_conditions = testCase.TestData.previous_conditions;
    
    %initial conditions
    Xn = [rand(M-1, 1)/1e+2; rand(M-1, 1)/1e+2; rand(M+1, 1)/1e+2; 1; 1];

    L = length(Xn);
    ftm = @(Xn) function_to_minimize_v1(Xn, previous_conditions, dt, contact_points, testCase.TestData.settings);
    expected_value = zeros(length(ftm(Xn)), length(Xn));
    COL = ftm(Xn);
    ddt = 1e-5;
    for jj = 1:L
        perturb = zeros(L, 1); perturb(jj) = ddt;
        expected_value(:, jj) = (ftm(Xn + perturb) - COL)/ddt;
    end
    actual_value = JacobianCalculator_v1(Xn, previous_conditions, dt, contact_points, testCase.TestData.settings);
    
    verifyEqual(testCase, actual_value, expected_value, 'AbsTol', 1e-5); 
end

function testJacobianContactPoints(testCase)
    % No contact points, zero initial conditions, zero testing value
    dt = 1;
    contact_points = 5;
    M = testCase.TestData.M;
    
    %current_conditions  = testCase.TestData.current_conditions;
    previous_conditions = testCase.TestData.previous_conditions;
    
    %initial conditions
    Xn = [rand(M-1, 1)/1e+2; rand(M-1, 1)/1e+2; rand(M+1, 1)/1e+2; 1; 1];

    L = length(Xn);
    ftm = @(Xn) function_to_minimize_v1(Xn, previous_conditions, dt, contact_points, testCase.TestData.settings);
    expected_value = zeros(length(ftm(Xn)), length(Xn));
    COL = ftm(Xn);
    ddt = 1e-5;
    for jj = 1:L
        perturb = zeros(L, 1); perturb(jj) = ddt;
        expected_value(:, jj) = (ftm(Xn + perturb) - COL)/ddt;
    end
    actual_value = JacobianCalculator_v1(Xn, previous_conditions, dt, contact_points, testCase.TestData.settings);
    
    verifyEqual(testCase, actual_value, expected_value, 'AbsTol', 1e-5); 
end

function testJacobianRandomContactPoints(testCase)
    % No contact points, zero initial conditions, zero testing value
    dt = 1;
    
    M = testCase.TestData.M;
    contact_points = floor(testCase.TestData.N/2 * rand());
    
    %current_conditions  = testCase.TestData.current_conditions;
    previous_conditions = testCase.TestData.previous_conditions;
    
    %initial conditions
    Xn = [rand(M-1, 1)/1e+2; rand(M-1, 1)/1e+2; rand(M+1, 1)/1e+2; rand()/1e+2; rand()/1e+2];

    L = length(Xn);
    ftm = @(Xn) function_to_minimize_v1(Xn, previous_conditions, dt, contact_points, testCase.TestData.settings);
    expected_value = zeros(length(ftm(Xn)), length(Xn));
    COL = ftm(Xn);
    ddt = 1e-5;
    for jj = 1:L
        perturb = zeros(L, 1); perturb(jj) = ddt;
        expected_value(:, jj) = (ftm(Xn + perturb) - COL)/ddt;
    end
    actual_value = JacobianCalculator_v1(Xn, previous_conditions, dt, contact_points, testCase.TestData.settings);
    
    verifyEqual(testCase, actual_value, expected_value, 'AbsTol', 1e-5); 
end

function testJacobianRandomPreviousConditions(testCase)
    % No contact points, zero initial conditions, zero testing value
    dt = 1;
    
    M = testCase.TestData.M;
    contact_points = floor(testCase.TestData.N/2.2 * rand());
    
    %current_conditions  = testCase.TestData.current_conditions;
    previous_conditions = testCase.TestData.previous_conditions;
    
    previous_conditions{1}.deformation_amplitudes = [0; rand(M-1, 1)/100];
    previous_conditions{1}.deformation_velocities = [0; rand(M-1, 1)/100]';
    previous_conditions{1}.pressure_amplitudes    = rand(M+1, 1)/100;
    previous_conditions{1}.center_of_mass         = rand()/100;
    previous_conditions{1}.center_of_mass_velocity= rand();

    previous_conditions{2}.deformation_amplitudes = [0; rand(M-1, 1)/100];
    previous_conditions{2}.deformation_velocities = [0; rand(M-1, 1)/100]';
    previous_conditions{2}.pressure_amplitudes    = rand(M+1, 1)/100;
    previous_conditions{2}.center_of_mass         = rand()/100;
    previous_conditions{2}.center_of_mass_velocity= rand();

    %initial conditions
    Xn = [rand(M-1, 1)/1e+4; rand(M-1, 1)/1e+4; rand(M+1, 1)/1e+4; rand()/1e+2; rand()/1e+2];

    L = length(Xn);
    ftm = @(Xn) function_to_minimize_v1(Xn, previous_conditions, dt, contact_points, testCase.TestData.settings);
    expected_value = zeros(length(ftm(Xn)), length(Xn));
    COL = ftm(Xn);
    ddt = 1e-5;
    for jj = 1:L
        perturb = zeros(L, 1); perturb(jj) = ddt;
        expected_value(:, jj) = (ftm(Xn + perturb) - COL)/ddt;
    end
    actual_value = JacobianCalculator_v1(Xn, previous_conditions, dt, contact_points, testCase.TestData.settings);
    
    verifyEqual(testCase, actual_value, expected_value, 'AbsTol', 1e-5); 
end


