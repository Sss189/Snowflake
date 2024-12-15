function snowflake_dla()
    % Parameters
    L = 401;                        % Grid size (odd for central seed)
    max_radius = 380;              % Maximum radius for growth
    center = (L + 1) / 2;          % Center of the grid
    num_particles = 600;           % Number of particles released
    grid = zeros(L);               % Initialize grid
    grid(center, center) = 1;      % Seed at the center
    axes_angles = (0:5) * pi / 3;  % Symmetry axes angles

    % Initialize maximum radius for each axis
    max_radius_axes = zeros(1, 6); 

    % Create figure
    figure;
    hold on;
    axis equal;
    axis off;

    % Simulation loop
    for n = 1:num_particles
        % Release particle from random position
        angle = 2 * pi * rand;
        x = round(center + max_radius * cos(angle));
        y = round(center + max_radius * sin(angle));

        is_attached = false;

        while ~is_attached
            % Random walk
            direction = randi(4);
            switch direction
                case 1, x = x + 1; % Move right
                case 2, x = x - 1; % Move left
                case 3, y = y + 1; % Move up
                case 4, y = y - 1; % Move down
            end

            % Check boundary conditions
            if sqrt((x - center)^2 + (y - center)^2) > max_radius
                angle = 2 * pi * rand;
                x = round(center + max_radius * cos(angle));
                y = round(center + max_radius * sin(angle));
            end

            % Check for attachment to growth site
            if x > 1 && x < L && y > 1 && y < L && ...
               (grid(x + 1, y) + grid(x - 1, y) + ...
                grid(x, y + 1) + grid(x, y - 1)) > 0
                
                % Calculate angle from center to particle
                dx = x - center;
                dy = y - center;
                theta = atan2(dy, dx);
                if theta < 0, theta = theta + 2 * pi; end
                
                % Calculate attachment probability
                angle_diff = abs(theta - axes_angles);
                [~, min_index] = min(angle_diff);
                attachment_probability = max(0.13, 1 - min(angle_diff) / (pi / 8));
                weight_factor = 3; % Attraction factor for main axis
                if angle_diff(min_index) < pi / 50
                    attachment_probability = attachment_probability * weight_factor;
                end

                if rand < attachment_probability
                    grid(x, y) = 1; % Attach particle
                    is_attached = true;
                    
                    % Update maximum radius for corresponding axis
                    current_radius = sqrt(dx^2 + dy^2);
                    max_radius_axes(min_index) = max(max_radius_axes(min_index), current_radius);
                    
                    % Add symmetric copies
                    grid = plot_symmetric(x, y, center, grid);
                    drawnow limitrate;
                end
                
                % Force attachment on main axes
                for m = 1:6
                    angle = axes_angles(m);
                    for r = 1:max_radius_axes(m)
                        % Calculate point on axis
                        x_force = round(center + r * cos(angle));
                        y_force = round(center + r * sin(angle));
                        if x_force > 1 && x_force < L && y_force > 1 && y_force < L && grid(x_force, y_force) == 0
                            grid(x_force, y_force) = 1; % Force attachment
                            grid = plot_symmetric(x_force, y_force, center, grid);
                        end
                    end
                end
            end
        end
    end
    
    % Generate filled grid and create 3D visualization
    create3DView(grid);
end

function grid = plot_symmetric(x, y, center, grid)
    axes_angles = (0:5) * pi / 3; % Symmetry axes
    [dx, dy] = deal(x - center, y - center);
    theta = atan2(dy, dx);
    if theta < 0, theta = theta + 2 * pi; end

    % Find nearest symmetry axis
    [~, idx] = min(abs(theta - axes_angles));
    alpha = axes_angles(idx);

    % Reflect across nearest symmetry axis
    cos2a = cos(2 * alpha);
    sin2a = sin(2 * alpha);
    dx_reflected = cos2a * dx + sin2a * dy;
    dy_reflected = sin2a * dx - cos2a * dy;

    hold on;
    for i = 0:5
        angle = i * pi / 3;
        cos_angle = cos(angle);
        sin_angle = sin(angle);
        
        % Rotate original point
        xn = center + cos_angle * dx - sin_angle * dy;
        yn = center + sin_angle * dx + cos_angle * dy;
        rectangle('Position', [round(xn) - 0.5, round(yn) - 0.5, 1, 1], ...
                  'FaceColor', [0.53, 0.81, 0.98], ...
                  'EdgeColor', 'none');
        grid(round(xn), round(yn)) = 1; 

        % Rotate reflected point
        xn_reflected = center + cos_angle * dx_reflected - sin_angle * dy_reflected;
        yn_reflected = center + sin_angle * dx_reflected + cos_angle * dy_reflected;
        rectangle('Position', [round(xn_reflected) - 0.5, round(yn_reflected) - 0.5, 1, 1], ...
                  'FaceColor', [0.53, 0.81, 0.98], ...
                  'EdgeColor', 'none');
        grid(round(xn_reflected), round(yn_reflected)) = 1; 
    end
end

function create3DView(grid)
    % Fill holes and apply Gaussian filter
    filled_grid = imclose(grid, strel('square', 2)); 
    filled_grid = imfill(filled_grid, 'holes'); 
    sigma = 0.4;
    filled_grid = imgaussfilt(double(filled_grid), sigma) > 0; 
    filled_grid = imfill(filled_grid, 'holes'); 

    % Create Z data for height
    Z = double(filled_grid) * 5; 
    [x, y] = meshgrid(1:size(filled_grid, 2), 1:size(filled_grid, 1));

    % Increase resolution
    factor = 4;
    [xq, yq] = meshgrid(linspace(1, size(filled_grid, 2), size(filled_grid, 2) * factor), ...
                        linspace(1, size(filled_grid, 1), size(filled_grid, 1) * factor));
    
    % Bilinear interpolation
    Zq = interp2(x, y, Z, xq, yq, 'linear');

    % Plot 3D snowflake structure
    figure('Units', 'Inches', 'Position', [2.8, 0.6, 7, 6]);
    surf(xq, yq, Zq, 'EdgeColor', 'none'); % No edge color
    colormap(cool); % Color map
    view(90, 89); 
    axis equal; 
    axis off;
    title('3D Snowflake Structure'); 
    
    % Enhance lighting
    light; 
    lighting phong; 
end