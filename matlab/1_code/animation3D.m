% Assuming r_data is an n_theta x n_t matrix with radial data for each theta at each time t
% theta_data is a 1D array with angular positions (n_theta)
% t_data is a 1D array with time values (n_t)

% Example Data Setup (use your actual data)
% r_data = ...; % your radius data
% theta_data = ...; % your angular positions in radians
% t_data = ...; % your time values

% Convert spherical (r, theta) to Cartesian coordinates (x, z) for plotting
[x_data, z_data] = pol2cart(theta_data, r_data); 

% Create the figure for animation
figure;
ax = axes('XLim', [-max(r_data(:)), max(r_data(:))], 'YLim', [-max(r_data(:)), max(r_data(:))], 'ZLim', [0, max(r_data(:))]);
hold on;
grid on;
view(3); % Set 3D view
xlabel('X');
ylabel('Y');
zlabel('Z');

% Droplet surface plot (initial state)
droplet_surface = surf(x_data, zeros(size(x_data)), z_data, 'EdgeColor', 'none');
colormap('turbo');
colorbar;

% Animation
for t = 1:length(t_data)
    % Update z data for the droplet
    z = r_data(:, t); % Update for current time step
    [x, z_cart] = pol2cart(theta_data, z); % Convert to Cartesian coordinates
    
    % Update surface plot
    set(droplet_surface, 'XData', x, 'YData', zeros(size(x)), 'ZData', z_cart);
    drawnow;
    
    % Optionally save frames for creating a video
    % frame = getframe(gcf);
    % writeVideo(your_video_writer, frame);
end
