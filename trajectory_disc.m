% *********************************************************************** %
%                                                                         %
% Project           : Golf Trajectory Simulator                           %
%                                                                         %
% File name         : trajectory_disc.m                                   %
%                                                                         %
% Version           : 1.0                                                 %
%                                                                         %
% Author            : Gabriele J. Giuli                                   %
%                                                                         %
% Date created      : 27/07/2021                                          %
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
    
    % Initial conditions
    xi = 0.8;
    yi = 0.8;
    dxi = 5;
    dyi = 0.7;
    
    % Computation parameters
    h = 1E-2;   % Step size
    t_f = 1;    % Final time
    
    % Surface curve generating function.
    f = @(x, y) x.^2 + y.^2;
    
    % Surface grid parameters
    x_i = -1;
    x_f = +1;
    x_h = 0.01;
    
    y_i = -1;
    y_f = +1;
    y_h = 0.01;
    
% -----------------------------------------------------------------------

% Generate surface mesh
x = x_i:x_h:x_f;
y = y_i:y_h:y_f;

% Calculate surface z-values
% This step will be replaced by the lidar data
[x_mtx, y_mtx] = meshgrid(x, y);
z_mtx = f(x_mtx, y_mtx);

% Find numerical gradient
[dfdx, dfdy] = gradient(z_mtx, x_h, y_h);

% Create RK function handles with linear interpolation
f1 = @(t, u1, u2, u3, u4) u2;
f2 = @(t, u1, u2, u3, u4) -g/m*interp2(x_mtx, y_mtx, dfdx, u1, u3) - a/m*u2;
f3 = @(t, u1, u2, u3, u4) u4;
f4 = @(t, u1, u2, u3, u4) -g/m*interp2(x_mtx, y_mtx, dfdy, u1, u3) - a/m*u4;

% Solve differential equations
[t_de, x_de, ~, y_de, ~] = RK4_4(f1, f2, f3, f4, 0, xi, dxi, yi, dyi, h, t_f);

% Set up vectors and scalars
N = length(t_de);
z_de = zeros(1, N);

% Find trajectory
for i = 1:1:N
    z_de(i) = interp2(x, y, z_mtx, x_de(i), y_de(i));
end

% Plot trajectory on surface
figure;
surf(x_mtx, y_mtx, z_mtx);
colormap summer
shading interp
hold on;
plot3(x_de, y_de, z_de, 'r', 'LineWidth', 3);
hold off;

% Create motion plot
figure;
plot(x_de, y_de, 'LineWidth', 1);
grid on;
box on;

% Create trajectory plot
figure;
plot3(t_de, x_de, y_de, 'LineWidth', 1);
grid on;
box on;

% Create time plots
figure;
subplot(2, 1, 1);
plot(t_de, x_de, 'b', 'LineWidth', 1);
subplot(2, 1, 2);
plot(t_de, y_de, 'r', 'LineWidth', 1)
