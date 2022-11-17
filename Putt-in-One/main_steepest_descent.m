% *********************************************************************** %
%                                                                         %
% Project           : Golf Trajectory Simulator                           %
%                                                                         %
% File name         : main_steepest_descent.m                             %
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
    %f = @(x, y) 0.1*cos(x);
    f = @(x, y) cos(x.^2 + y.^2) + 2*sin((x - 0.3).^2 + (y - 0.6).^2);
    
    % Scene geometry
    x_b = 0.9;  % Ball x-coordinate
    y_b = 0.9;  % Ball y-coordinate
    x_h = 0.2;  % Hole x-coordinate
    y_h = -0.3; % Hole y-coordinate
    r_h = 0.01; % Hole radius
    
    % Solving parameters
    h = 1E-2;   % Step size
    tau = 3;    % Number of time contants
    
    % Central difference stepsize
    cd_h = 0.001;
    gamma = 10;
    
% -----------------------------------------------------------------------

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

% Find numerical gradient and final execution time
[dfdx, dfdy] = gradient(z_mtx, x, y);
t_f = tau*m/a;

% Create initial guess.
% NOTE: The speed is corrected for the equations of motions solved
%       on a flat plane and to take into account the difference in 
%       potential energy for the ball at the starting and end position.
vel_abs = a/m*norm(bth_v) + sqrt(2*g*abs(z_h - z_b))*sign(z_h - z_b);
sol = vel_abs*bth_n;

% Set distance d > r_h for entry point in while loop
d = r_h + 1;

% Apply (local) gradient descent method to find local minimum
while (d > r_h)
    
    % Solve iteratively
    [sol, d, x_de, y_de] = iterate_solution(sol, cd_h, gamma, x_h, ...
        y_h, x_mtx, y_mtx, dfdx, dfdy, x_b, y_b, h, t_f, g, m, a);
     
    % Plot solution
    pcolor(x_mtx, y_mtx, z_mtx);
    colormap summer
    shading interp
    hold on;
    plot(x_b, y_b, 'w.', 'MarkerSize', 15);
    plot(x_h, y_h, 'k.', 'MarkerSize', 15);
    plot(x_de, y_de, 'r', 'LineWidth', 2);
    title(d);
    hold off;
    
    % Wait for plot to render
    pause(0.2);
        
end

% Create z-values
z_de = f(x_de, y_de); 

% Plot 3D scene
figure;
surf(x_mtx, y_mtx, z_mtx);
colormap summer
shading interp
hold on;
plot3(x_b, y_b, z_b, 'w.', 'MarkerSize', 15);
plot3(x_h, y_h, z_h, 'k.', 'MarkerSize', 15);
plot3(x_de, y_de, z_de, 'LineWidth', 3);
hold off;
