% *********************************************************************** %
%                                                                         %
% Project           : Golf Trajectory Simulator                           %
%                                                                         %
% File name         : main_trial_and_error.m                              %
%                                                                         %
% Version           : 1.0                                                 %
%                                                                         %
% Author            : Gabriele J. Giuli                                   %
%                                                                         %
% Date created      : 28/07/2021                                          %
%                                                                         %
% *********************************************************************** %

close all;
clear all;
clc;

% -------------------------> GLOBAL PARAMETERS <-------------------------
    
    % Physical constants
    g = 9.81;   % Gravitational acceleration
    m = 0.1;    % Golf ball mass
    a = 0.6;    % Friction coefficient
    
    % x-axis quantization parameters
    x_i = -1;   % x mimimum value
    x_f = +1;   % x maximum value
    x_s = 0.01; % x quantization
    
    % y-axis quantization parameters
    y_i = -1;   % y minimum value
    y_f = +1;   % y maximum value
    y_s = 0.01; % y quantization
    
    % Surface curve generating function.
    %f = @(x, y) x.^2 + y.^2 + x;
    %f = @(x, y) 0.1*x + 0.01*y;
    f = @(x, y) 0.1*cos(x);
    
    % Scene geometry
    x_b = 0.9;  % Ball x-coordinate
    y_b = 0.9;  % Ball y-coordinate
    x_h = 0.2;  % Hole x-coordinate
    y_h = -0.3; % Hole y-coordinate
    r_h = 0.05; % Hole radius
    
    % Simulation parameters
    delta_a = -10:1:10;
    delta_v = 10:1:12; 
    
    % Solving parameters
    h = 1E-2;   % Step size
    t_f = 1;    % Final time
    
% -----------------------------------------------------------------------

% Initialize waitbar
count_a = length(delta_a);
count_v = length(delta_v);
total = count_a*count_v;
current = 1;
wb = waitbar(0, 'Initializing...');

% Generate surface mesh
x = x_i:x_s:x_f;
y = y_i:y_s:y_f;

% Calculate surface z-values
% This step will be replaced by the lidar data
[x_mtx, y_mtx] = meshgrid(x, y);
z_mtx = f(x_mtx, y_mtx);

% IMPORTANT: Below this point it is no longer possible to call the
%            surface generating function f(x, y) in order to mimic
%            the final behaviour of the program (i.e. data loaded)
%            from LIDAR scan.

% Find ball and hole z-coordinates
z_b = interp2(x_mtx, y_mtx, z_mtx, x_b, y_b);
z_h = interp2(x_mtx, y_mtx, z_mtx, x_h, y_h);

% Find ball-to-hole normal vector
bth_v = [x_h - x_b; y_h - y_b];
bth_n = bth_v/norm(bth_v);

% Find numerical gradient
[dfdx, dfdy] = gradient(z_mtx, x, y);

% Iterate through angles and velocity magnitudes
for theta = delta_a*pi/180
    for v = delta_v
        
        % Update waitbar
        waitbar(current/total, wb, 'Processing trajectories...');
        current = current + 1;
        
        % Find rotation matrix for theta
        R = [cos(theta) -sin(theta);
             sin(theta) cos(theta)];
         
        % Calculate initial speed vector from rotation and scaling
        throw_v = v*R*bth_n;
        
        % Calculate trajectory
        % NOTE: the parameters x_mtx and y_mtx could be substituted
        %       with x and y respectively. This could prove useful
        %       later on in the development.
        [t_de, x_de, y_de] = solve_trajectory(x_mtx, y_mtx, dfdx, dfdy, ...
            x_b, throw_v(1), y_b, throw_v(2), h, t_f, g, m, a);
        
        % Check if trajectory goes in hole
        d = calculate_distance(x_de, y_de, x_h, y_h);
        
        % Check if ball goes in hole
        if d <= r_h
            theta
            v
        end
        
        % Plot figure
        pcolor(x_mtx, y_mtx, z_mtx);
        colormap summer
        shading interp
        hold on;
        plot(x_b, y_b, 'w.', 'MarkerSize', 15);
        plot(x_h, y_h, 'k.', 'MarkerSize', 15);
        plot(x_de, y_de, 'r', 'MarkerSize', 15);
        hold off;
        
    end
end

% Close waitbar
close(wb);

% Plot 3D scene
figure;
surf(x_mtx, y_mtx, z_mtx);
colormap summer
shading interp
hold on;
plot3(x_b, y_b, z_b, 'w.', 'MarkerSize', 15);
plot3(x_h, y_h, z_h, 'k.', 'MarkerSize', 15);
hold off;
