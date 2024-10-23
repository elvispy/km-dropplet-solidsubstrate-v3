function h = plot_condition(idx, conditions, varargin)
%{
    This function will plot the current conditions of the problem.

    INPUTS:
      - idx: The position index of the plot. Must be 1 or 2, so as to
      position the box plot.
      - conditions: current conditions to be plotted. Can be a
      ProblemConditions_vx (x<=3) struct. conditions also can be a
      one-dimensional scalar vector representing the legendre amplitudes,
      and contact will be assumed
    Optional inputs:
      - varargin{1} = N: scaling factor for plotting. If radius of the
      undisturbed sphere is 1, then it wil show from -N to N. Default is
      two
      - varargin{2}: Plot title
      - varargin{3}: Angles vector to determine contact angle

%}
    h = figure(idx);
    if idx == 2
        set(gcf, 'Position', [780 159 760 586]);
    else
        set(gcf, 'Position', [0   159 760 586]);
    end
    clf;
    hold on;  
    cut = 0.75 * pi;
    sample = [linspace(0, cut, 30), linspace(cut, pi, 30)];
    arrX = sin(sample);
    arrY = cos(sample);
    etas = zeta_generator(conditions);
    
    if isstruct(conditions)
        height = conditions.center_of_mass;
        %plot([-conditions.contact_radius, conditions.contact_radius], [0, 0], 'g--', 'LineWidth', 3);
        %yline(, 'y', 'LineWidth', 2);
        nb_harmonics = conditions.nb_harmonics;
    else
        nb_harmonics = length(conditions);
        height = 1 + sum(arrayfun(@(idx) (-1)^idx * conditions(idx), 2:nb_harmonics));
    end
    EtaX = arrayfun(@(angle) sin(angle) * (1+  etas(angle)), sample);
    EtaY = height + arrayfun(@(angle) cos(angle) .* (1+  etas(angle)), sample);
    if nargin >= 5
        theta_vector = varargin{3}.theta_vector;
        pressure_unit = varargin{3}.pressure_unit;
        EtaX2 = arrayfun(@(angle) sin(angle) * (1+  etas(angle)), theta_vector);
        EtaY2 = height + arrayfun(@(angle) cos(angle) .* (1+  etas(angle)), theta_vector);
        scatter([EtaX2, -EtaX2], [EtaY2, EtaY2], 30, 'filled');
        if conditions.contact_points > 0
            contact_angle = theta_vector(conditions.contact_points);
            etacontactX = sin(contact_angle) * (1 + etas(contact_angle));
            etacontactY = cos(contact_angle) * (1 + etas(contact_angle)) + height;
            
            %scatter(-EtaX2, EtaY2, 50, 'filled');
            scatter(etacontactX, etacontactY, 60, 'Marker', 'x', 'MarkerEdgeColor', 'r', 'LineWidth', 5);
            
            contact_radius = sin(theta_vector(conditions.contact_points)/2 + theta_vector(conditions.contact_points+1)/2) *...
                (1 + etas(theta_vector(conditions.contact_points)/2 + theta_vector(conditions.contact_points+1)/2));
            plot([-contact_radius, contact_radius], [0, 0], 'g--', 'LineWidth', 3);
        end
    else
        pressure_unit = 1;
        %warning("Assuming pressure unit");
    end
    
    fill( EtaX,EtaY, [135, 206, 235]/256, 'LineStyle','none' ,'FaceAlpha', 0.3);
    fill(-EtaX,EtaY, [135, 206, 235]/256, 'LineStyle','none', 'FaceAlpha', 0.3);
    
%     gr = -100:100;
%     if nargin > 2
%         dr = varargin{1}; 
%         scatter(dr * gr, 0 * gr, 'Marker', 'o', 'MarkerEdgeColor', 'r', 'LineWidth', 2, 'SizeData', 10 * ones(size(gr)));
%     end
    
    %scatter(0, 0, 'Marker', 'o', 'MarkerEdgeColor', 'b', 'LineWidth', 2, 'SizeData', 10);
    if isstruct(conditions)
        zps = @(theta) sum(conditions.pressure_amplitudes(2:end)' .* collectPl(length(conditions.pressure_amplitudes(2:end)), ...
            cos(theta)), 1);
        %zps = zeta_generator(conditions.pressure_amplitudes((end-nb_harmonics+1):end));
        ps = @(ang) zps(ang) + conditions.pressure_amplitudes(1);
        mps = ps(sample)/pressure_unit;

        %mps(5) = 1;
        quiver(EtaX, EtaY, mps .* (-arrX), mps .* (-arrY), 'AutoScaleFactor',1);

        if idx ~= 1
            angle_tol = pi*2/nb_harmonics;
            cm = conditions.center_of_mass;
            plot([0, sin(pi+angle_tol)], [cm, cos(pi+angle_tol)], 'b','LineWidth',0.3);
            plot([0, sin(pi-angle_tol)], [cm, cos(pi-angle_tol)], 'b','LineWidth',0.3);
        end
    end
    
    if nargin >= 5 && isstruct(conditions)
        theta_vector = varargin{3}.theta_vector;
        if conditions.contact_points > 0
            x = 0.15 + ceil(3 * conditions.center_of_mass * sin(pi - theta_vector(conditions.contact_points)))/3;
        else
            x = 0.55;
        end
        %x = 0.15 + ceil(3*conditions.contact_radius)/3; % 0.1 + conditions.contact_radius;
    else
        x = 0.55;
    end
    if nargin > 2
        N = varargin{1};
    else
        N = 2;
    end
    x = 1;
    xlim([-N*x, N*x]);
    ylim([-2*N/5 * x , 0.8*2*N*x ]);
    %axis equal;
    %ylim([-1.5, 1.5]);
    yline(0, 'k', 'LineWidth', 1.5);
    
    if nargin > 3 && idx == 2 && isstruct(conditions)
        ss = varargin{2};
        
        %s = sprintf("Attempting to fit the solution with contact radius %.3f. \n Previous contact radius: %.3f. Iteration number: %d", ...
        %    conditions.contact_radius, ss.contact_radius, ss.iteration);
        %title(s, 'FontSize', 14);
        legend("Contact Radius");
        x = xlim;
        y = ylim;
        y = y(2) - (y(2) - y(1))/10;
        x = x(1) + (x(2) - x(1))/10;
        text(x, y, sprintf("v_{cm} = %.7g", conditions.center_of_mass_velocity), 'FontSize', 14);
        text(x, y - y/10, sprintf("z_{cm} = %.7g", conditions.center_of_mass), 'FontSize', 14);
        %pause(0);
    elseif nargin >= 4
        
        title(varargin{2}, 'FontSize', 14);
        drawnow limitrate;
    end
    
end