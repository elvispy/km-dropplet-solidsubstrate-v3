
close all; clear all;
addpath(fullfile(pwd, "..", "..", "1_code", "simulation_code"));
curr = pwd;
%p = fullfile(pwd, "..", "D50Quant100", "rho1000sigma7220nu98muair0", "RhoS1000SigmaS7220", "R0350mm", "ImpDefCornerAng180U39", "N=20tol=5.00e-05"); % uigetdir();
%cd(p);

saving = false;

myFont = "Arial";
n=8000;
amps = [0, -0.144, 0.0695, -0.0175, 0.0138, -8.13e-04, -2.37e-03, 4.05e-03, -4.73e-03, 3.59e-03, -2.15e-03, 8.12e-04, 2.23e-04, -9.30e-04, 1.27e-03, -1.19e-03, 8.69e-04, -4.11e-04, -3.66e-05, 3.64e-04, -5.31e-04, 5.42e-04, -4.25e-04, 2.28e-04, -1.15e-05, -1.69e-04, 2.76e-04, -2.97e-04, 2.44e-04, -1.40e-04, 1.83e-05, 8.97e-05, -1.60e-04, 1.80e-04, -1.54e-04, 9.31e-05, -1.79e-05, -5.18e-05, 9.96e-05, -1.17e-04, 1.03e-04, -6.49e-05, 1.57e-05, 3.15e-05, -6.53e-05, 7.89e-05, -7.12e-05, 4.67e-05, -1.32e-05, -1.98e-05, 4.43e-05, -5.49e-05, 5.06e-05, -3.41e-05, 1.08e-05, 1.28e-05, -3.06e-05, 3.89e-05, -3.64e-05, 2.49e-05, -8.40e-06, -8.53e-06, 2.15e-05, -2.76e-05, 2.60e-05, -1.79e-05, 6.12e-06, 5.98e-06, -1.52e-05, 1.95e-05, -1.82e-05, 1.23e-05, -3.87e-06, -4.63e-06];

thetas=0:2*pi/n:2*pi;
f = zeta_generator(amps);
r = 1 + f(thetas); 
x = r .* sin(thetas);
CENTER = 0.7686;
y = r .* cos(thetas);
saving_figure = figure;
plot(x,y + CENTER,'Color', [.1 .08 .08] ,'LineWidth',3);
myax = gca; mypos = myax.Position;
axis equal
hold on
mylims = 1.75;
set(myax,'Xlim',[-mylims mylims],'Xtick',[],'Ytick',[]);
set(myax, 'Ylim', [-.4, -1.4+2*mylims]);

%S(t)
epsi = 1;
idx1 = floor((4-epsi)/8*n):floor((4+epsi)/8*n);
x1 = x(idx1);
y1 = y(idx1);
%plot(x1,y1-.02,'color',[.5, .5, .5], 'LineWidth', 2);%[0/256, 191/255, 255/255],'LineWidth',4)

%[xSt, ySt] = gca_to_Normalized(myax, [x1(end) + 0.02, x1(end) + 0.29], ...
%    [y1(1) - 0.28, y1(1)-0.30]);
%annotation('textarrow', xSt, ySt, 'HeadLength', 10, ...
%'HeadStyle', 'vback3', 'String', '$S(t)$', 'Interpreter', 'latex', ...
%'LineWidth', 1.5, 'FontSize', 14);
hold on


% Free surfaces (eta)
t2 = (-0.01):.01:mylims;
x2 = t2; %t2+x1(length(x1));
%y2 = y1(length(y1))+(1-exp(-t2/2.5))-.02;
ZERO = 0;

A = 0.5; k = 3; alpha = 1.5; omega = 0.5;
r = linspace(0, pi, length(x2));  % Radial distance
z =  -A * sin(k * r) .* (1-exp(-alpha * (r-0.5).^2)) .* (1-exp(-0.1 * r.^2));

%cc = plot(x2,z,'color',[.5, .5, .5],'LineWidth',2);%free surfaces

%x3 = -x2;
%plot(-x2,z,'color',[.5, .5, .5],'LineWidth',2)%free surfaces
hold on

% z = 0
%plot([-.5 .5], [ZERO, ZERO], '--','color',[.8 .0 .08],'LineWidth',4);%free surfaces
aa = plot([-.53 .53], [ZERO, ZERO], '--','color',[.8 .0 .08],'LineWidth',4);%free surfaces
cc = fill([-10, 10, 10, -10], [-0.0005, -0.0005, -10, -10], [235 176 0]./256, 'EdgeColor', 'none', 'FaceAlpha',0.3);
% h(t)
%Centre mark

%plot([-val, val], [val+disloc, -val+disloc], 'k');
%plot([-val, val], [-val+disloc, val+disloc], 'k');

%plot([-0.02 -.15],[disloc disloc], 'Color', [.7 .7 .7], 'LineWidth', 1.1);%bottom measure transporter for h
[xht, yht] = gca_to_Normalized(myax, [sin(3*pi/4)*.99+.15, ZERO], [(cos(3*pi/4)+1)*.99 + CENTER+.15, ZERO+CENTER]);
annotation('doublearrow',xht, yht, ...
    'Head1Length', 9, 'Head2Length', 0, ...
    'Head1Style', 'vback3', 'Head2Style', 'vback3' ,'Linewidth',2.5 )
text(sin(3*pi/4)*0.33, CENTER, "$\zeta(\theta, t)$", ...
    'Fontsize', 20, 'FontName', myFont, 'interpreter', 'latex')

% Plot angle arc using annotation for the angle (phi)
theta = linspace(0, pi/3.1, 20);
x_arc = 0.10*sin(theta);
y_arc = 0.10*cos(theta) + CENTER;
plot(x_arc, y_arc, 'k', 'LineWidth', 1.5);

% Text annotation for the angle \varphi
text(ZERO +.05, ZERO + CENTER + 0.25, '$\theta$', 'FontSize', 20, 'Interpreter', 'latex');
%annotation('doublearrow', x_arc(end-1:end), y_arc(end-1:end), 'Head1Length', 9, 'Head2Length', 0, ...
%  'Head1Style', 'vback3', 'Linewidth', .75);
[xht, yht] = gca_to_Normalized(myax, [sin(0)*.99, ZERO], [(cos(0))*.99 + CENTER - 0.15, ZERO+CENTER]);
annotation('doublearrow',xht, yht, ...
    'Head1Length', 0, 'Head2Length', 0, ...
    'Head1Style', 'vback3', 'Head2Style', 'vback3' ,'Linewidth',2.5 );


% FLuid properties
%text(-1.34, -0.2, "$(\sigma, \rho, \nu)$", ...
%    'Fontsize', 20, 'FontName', myFont, 'interpreter', 'latex')
text(-0.84, 0.9, "$(\sigma, \rho, \nu)$", ...
    'Fontsize', 20, 'FontName', myFont, 'interpreter', 'latex')



legend([cc, aa], ["Substrate", "Contact area"], ...
    'FontSize', 16,  'interpreter', 'latex', 'AutoUpdate', 'off');
%legend([Ct, Lt, aa], ["$C(t)$", "$L(t)$", "$z = 0$"], ...
%    'FontSize', 11,  'interpreter', 'latex', 'AutoUpdate', 'off');
cd(curr);

if saving == true
    saveas(saving_figure, "schematicsdfd.png");
    %print(saving_figure, '-depsc', '-r300', "schematicsdfd.eps");
end

function [X, Y] = gca_to_Normalized(ca, xs, ys)
    % (xs(i), ys(i)) represents a point. For annotation function, xs(1),
    % ys(1) represent the initial point of the annotation and xs(2), ys(2)
    % the final point. 
    pos = ca.Position;
    xLims = ca.XLim;
    yLims = ca.YLim;
    X = zeros(size(xs));
    Y = zeros(size(ys));
    for i = 1:length(xs)
       pctX = (xs(i) - xLims(1))/(xLims(2) - xLims(1)); 
       pctY = (ys(i) - yLims(1))/(yLims(2) - yLims(1));
       
       X(i) = pos(1) + pctX * pos(3);
       Y(i) = pos(2) + pctY * pos(4);
    end
end
% 
% function load_vars(str)
%     global errored
% 
%     if errored == true
%         str = "errored_" + str; 
%     end
% 
%     vars = load(str);
%     fn = fieldnames(vars);
%     for ii = 1:length(fn)
%         assignin('caller', fn{ii}, vars.(fn{ii}));
%     end
% 
% end
