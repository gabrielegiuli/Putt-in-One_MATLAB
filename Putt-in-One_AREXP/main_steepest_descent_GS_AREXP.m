% *********************************************************************** %
%                                                                         %
% Project           : Golf Trajectory Simulator                           %
%                                                                         %
% File name         : main_steepest_descent_GS.m                          %
%                                                                         %
% Version           : 1.0                                                 %
%                                                                         %
% Author            : Gabriele J. Giuli                                   %
%                                                                         %
% Date created      : 29/07/2021                                          %
%                                                                         %
% *********************************************************************** %

close all;
clear all;
clc;

% -------------------------> GLOBAL PARAMETERS <-------------------------

    % Data parameters
    file_name = "fake_green.ply";   % Data filename
    y_step = 0.05;                   % x-variable quantization
    x_step = 0.05;                   % y-variable quantization
    
    % Physical constants
    g = 9.81;   % Gravitational acceleration
    m = 0.1;    % Golf ball mass
    a = 0.001;  % Friction coefficient
    
    % Scene geometry
    x_b = -2;  % Ball x-coordinate
    y_b = 2;  % Ball y-coordinate
    x_h = 0.2;  % Hole x-coordinate
    y_h = -4; % Hole y-coordinate
    r_h = 0.01; % Hole radius
    
    % Solving parameters
    h = 1E-2;   % Step size
    tau = 3;    % Number of time contants
    
    % Central difference stepsize
    cd_h = 0.001;
    gamma_0 = 10;
    gamma_f = 0.1;
    tau_gs = 5;
    
% -----------------------------------------------------------------------

% Read data and display it
ptCloud = pcread(file_name);

% Load data into vectors for interpolation
x_int = double(ptCloud.Location(:, 1));
y_int = double(ptCloud.Location(:, 3));
z_int = double(ptCloud.Location(:, 2));

% Perform interpolation
% NOTE: It might be better to use TriScatteredInterp instead. Need
%       to do some more research in this topic
surf_func = scatteredInterpolant(x_int, y_int, z_int);

% Create Mesh
x = double(ptCloud.XLimits(1)):x_step:double(ptCloud.XLimits(2));
y = double(ptCloud.ZLimits(1)):x_step:double(ptCloud.ZLimits(2));
[x_mtx, y_mtx] = meshgrid(x, y);

% Fill in data
z_mtx = surf_func(x_mtx, y_mtx);

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
n = 0;

% Apply (local) gradient descent method to find local minimum
while (d > r_h)
    
    % Find gamma
    gamma = (gamma_0 - gamma_f)*exp(-n/tau_gs) + gamma_f;
    n = n + 1;
    
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
z_de = surf_func(x_de, y_de); 

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
