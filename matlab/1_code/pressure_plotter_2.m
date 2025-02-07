function pressure_plotter_2(varargin)
    % Add functions to calculate maximum width
    close all;
    safe_folder = fullfile(fileparts(mfilename('fullpath')), "simulation_code");
    addpath(safe_folder, '-begin');
    [file, path] = uigetfile("*.mat");
    fullfilepath = fullfile(path, file);
    %clear is_adim
    load(fullfilepath, "recorded_conditions", "default_physical", "PROBLEM_CONSTANTS");
    load(fullfilepath, "time_unit", "length_unit");
    %MM = 200; 
    %H = floor(linspace(1, recorded_conditions{1}.nb_harmonics, 20));
    %pressure_amplitudes = zeros(length(H), MM);
    %deformation_amplitudes = zeros(length(H), MM);
    %times = zeros(MM, 1);
    % Defining new pressure unit
    theta_vector = PROBLEM_CONSTANTS.theta_vector;
    pressure_unit = default_physical.rhoS * default_physical.initial_velocity.^2;
    nb_harmonics = recorded_conditions{1}.nb_harmonics;
    
    pressure_amplitudes = cell2mat(cellfun(@(x) x.pressure_amplitudes, ...
        recorded_conditions, 'UniformOutput',false))'/pressure_unit;
    deformation_amplitudes = cell2mat(cellfun(@(x) x.deformation_amplitudes, ...
        recorded_conditions, 'UniformOutput',false))'/length_unit;
    times = cellfun(@(x) x.current_time, recorded_conditions)'/time_unit;
    numl = cellfun(@(x) x.contact_points, recorded_conditions)';

    dr = 1/400;
    N = floor(length(times)*0.95); M = floor(1.5/dr);
    pfield_radial = zeros(M, N); cRadii = zeros(N, 1);
    indexes = floor(linspace(1, length(times), N));
    for ii = 1:N
        jj = indexes(ii);
        f = zeta_generator(pressure_amplitudes(:, jj));
        cRadius = r_from_spherical(theta_vector(max(numl(jj), 1)), deformation_amplitudes(:, jj));
        nbPoints = 0:dr:cRadius;
        thetaplots = theta_from_cylindrical(0:dr:cRadius, deformation_amplitudes(:, jj));
        pfield_now = f(thetaplots); %- sum(pressure_amplitudes(:, jj));
        pmean = mean(f(linspace(0, pi/10, 20))); % - sum(pressure_amplitudes(:, jj));
        
        pfield_radial(1:length(nbPoints), ii) = -pfield_now'+pmean;
        cRadii(ii) = cRadius;
    end

    % Create a figure

    % Plot the matrix using pcolor
    NN = ceil(1/dr);
    
    disp(sum(pfield_radial > 5, "all")/sum(pfield_radial<inf, "all")); 
    %disp(sum(pfield_radial < 0, "all")/sum(pfield_radial<inf, "all"));
    pfield_radial(pfield_radial > 5) = 5;
    pfield_radial(pfield_radial <= 0) = 0;
    p = pcolor(dr * (0:(M-1)), times(indexes) , pfield_radial');
    hold on
    plot(cRadii, times(indexes), 'r--', 'LineWidth',1)
    % Enhance the plot
    p.EdgeColor = 'none';            % Remove gridlines for a smoother look
    jet2 = jet; jet2(1, :) = 1;
    colormap(flipud(bone));                   % Choose a visually striking colormap (jet, parula, etc.)
    cb = colorbar;                        % Add a colorbar for reference
    shading interp;                  % Interpolate shading for a smoother gradient
    ylabel(cb, '$\frac{p}{\rho V_o^2}$','interpreter','Latex','FontName','Times',...
        'FontSize',20)
    % Add labels and title
    set(gca,'FontName','Times','FontSize',14);
    xlabel('   $x/R_o$   ','interpreter','Latex','FontName','Times','FontSize',20)
    ylabel('$t/t_{\sigma}$','interpreter','Latex','FontName','Times',...
        'FontSize',20,'rotation',90)
    %xlim([0, 1]); ylim([0, 5]);

    return
    saving_figure_pressure = figure('Position', [100, 50, 1000, 300]); % Wider than tall
    hold on;
    %recorded_conditions{ii}.amplitude_defor = 1;
    idxs = 1:(nb_harmonics+1);
    idxs = idxs(max(abs(pressure_amplitudes), [] , 2) < 700); %2.^(1:floor(log2(nb_harmonics)));
    %cmap = colormap('spring'); %disp(size(cmap));
    %numColors = size(cmap, 1);
    lol = jet(length(idxs)+2);
    colororder(flipud(lol(2:end, :)));%cmap(floor(linspace(1, numColors, length(idxs))), :));
    %disp(size(times)); disp(size(deformation_amplitudes));
    plot(times, pressure_amplitudes(idxs(1:(end-1)), :), 'LineWidth',2);
    % Only show a subset of 
    set(gca, 'FontSize', 16);
    legend(arrayfun(@(i) string(i), idxs), 'FontSize', 12);
    xlabel('$ t/t_s $', 'Interpreter','Latex', 'FontSize', 20);
    ylabel('$ \frac{p}{\rho V_o^2} $', 'Interpreter','Latex', ...
        'Rotation', 90,  'FontSize', 20);
    sigmaS = default_physical.sigmaS; rhoS = default_physical.rhoS; Ro = default_physical.undisturbed_radius;
    Westar = default_physical.rhoS * default_physical.initial_velocity^2 * ...
        default_physical.undisturbed_radius / default_physical.sigmaS;
    Oh = default_physical.nu * sqrt(default_physical.rhoS / (default_physical.sigmaS ...
        * default_physical.undisturbed_radius));
    Bo = rhoS * default_physical.g * Ro^2 / sigmaS;

    title(sprintf("Deformation Amplitudes with We = %.3g, Oh = %.3g, Bo = %.3g", Westar, Oh, Bo));
    %yl = get(gca, 'YLim'); yl = [-max(abs(yl)), max(abs(yl))];
    %set(gca, 'YLim', yl); set(gca, 'XLim', 1.05*get(gca, 'XLim'))
    
end



%title('Contact radius and pressure field evolution');

% Adjust axes
%axis tight;                      % Fit the plot closely around the data

% cd(curr);
% saveas(saving_figure, "../../0_data/manual/pressure_plotter", 'fig');
% print(saving_figure, '-depsc', '-r300', "../../0_data/manual/pressure_plotter.eps");


function r = r_from_spherical(angles, amplitudes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r_from_spherical.m - Returns cylindrical coordinates from spherical coordinates
%
% Calculate the radius in cylindrical coordinates that
% corresponds to the azimutal angle in spherical coordinates
% given certain amplitudes. (angle pi points downwards)
% r and angle are related by the relation
% r = sin(angle * (1 + \sum_{i=1}^{N}  A_l(i) * Pl(cos(angle)))
% where Pl is the l-th legendre Polynomial. (see
% https://en.wikipedia.org/wiki/Legendre_polynomials#Rodrigues'_formula_and_other_explicit_formulas) 

% Arguments:
% - angles: One dimensional scalar array such that size(r) = size(angles)
% angles(i) corresponds to the angle that gives radius r(i).
% - A_l: One dimensional scalar array or STRUCT with field
% "deformation_amplitudes"
%
% Outputs:
% - r: One dimensional scalar array of non-negative entries
%
% EXAMPLES
% r_from_spherical(pi, rand(10, 1));   % Returns 0.0
%
% Written by: Elvis Aguero- 01/01/2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isstruct(amplitudes); amplitudes = amplitudes.deformation_amplitudes; end
    if size(amplitudes, 2) > 1; amplitudes = amplitudes'; end
    r = sin(angles) .* (1 + sum(amplitudes .* collectPl(length(amplitudes), cos(angles)), 1));
%     r = arrayfun(@(angle) ...
%         sin(angle) * (1 + sum(dot(amplitudes, ...
%         arrayfun(@(idx) LP{idx}(cos(angle)), 1:d)))), theta);
end

function angle = theta_from_cylindrical(r, A_l)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theta_from_cylindrical.m - Returns spherical coordinates from cylindrical coordinates
%
%  Calculate the azimutal angle in spherical coordinates
% that corresponds to radius r in cylindrical coordinates, 
% given amplitudes A_l. (angle pi points downwards)
% This algorithm uses Newton-Raphson to approximate the angle
% r and angle are related by the relation
% r = sin(angle) * (1 + \sum_{i=1}^{N}  A_l(i) * Pl(cos(angle)))
% where Pl is the l-th legendre Polynomial. (see
% https://en.wikipedia.org/wiki/Legendre_polynomials#Rodrigues'_formula_and_other_explicit_formulas) 

% Arguments:
% - r: One dimensional scalar array of non-negative entries
% - A_l: One dimensional scalar matrix or STRUCT with field
% "deformation_amplitudes"
%
% Outputs:
% - angles: One dimensional scalar array such that size(r) = size(angles)
% angles(i) corresponds to the angle that gives radius r(i).

% EXAMPLES
% theta_from_cylindrical(0, rand(10, 1));   % Returns pi
%
% Written by: Elvis Aguero- 01/01/2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isstruct(A_l); A_l = A_l.deformation_amplitudes; end
    if size(A_l, 2) > 1; A_l = A_l'; end

    zeta = zeta_generator(A_l);
    
    % Derivative of the function
    f_prime = @(theta) cos(theta) .* (1 + zeta(theta)) - sin(theta).^2 .* sum(times(A_l, collectdnPl(length(A_l), cos(theta))), 1);
    
    
    angle = zeros(size(r));

    % Newton Method!
    for ii = 1:length(r)
        % Function to be minimized
        if r(ii) == 0; angle(ii) = pi; continue; end
        f_objective = @(theta) sin(theta) .* (1 + zeta(theta)) - r(ii);

        if ii>=5 
            theta = interp1(r((ii-4):(ii-1)), angle((ii-4):(ii-1)), r(ii), 'makima', 'extrap');
        else
            theta = pi - asin(min(1, r(ii)));
        end
        tol_theta = 1e-7;
        n = 1;
        
        while abs(f_objective(theta)) >= tol_theta && n < 350
            theta = mod(theta - f_objective(theta)/f_prime(theta) - 1e-4, pi/2) + 1e-4 + pi/2; % If solution is close to pi, theta is unstable with mod function (therefore 1e-4 added)
            n = n + 1;
            if n == 200
                theta = 3.141592653589793;
            elseif n == 300
                theta = pi - asin(min(1, r(ii)));
                %theta = rand() * pi/2 + pi/2;
            end
        end
        if n >= 340; warning("Theta_from_cylindrical didnt converge"); end
        angle(ii) = theta;
    end

end