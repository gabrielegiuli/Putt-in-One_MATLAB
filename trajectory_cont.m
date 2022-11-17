% *********************************************************************** %
%                                                                         %
% Project           : Golf Trajectory Simulator                           %
%                                                                         %
% File name         : trajectory_cont.m                                   %
%                                                                         %
% Version           : 1.0                                                 %
%                                                                         %
% Author            : Gabriele J. Giuli                                   %
%                                                                         %
% Date created      : 19/07/2021                                          %
%                                                                         %
% *********************************************************************** %

close all;
clear all;
clc;

% Define variables for surface 
syms x y u2n u4n;

% -------------------------> GLOBAL PARAMETERS <-------------------------
    
    % Physical constants
    g = 9.81;   % Gravitational acceleration
    m = 0.1;    % Golf ball mass
    a = 0.6;    % Friction coefficient
    
    % Initial conditions
    xi = 0.1;
    yi = 0.1;
    dxi = 0.1;
    dyi = 0.7;
    
    % Computation parameters
    h = 1E-3;   % Step size
    t_f = 5;    % Final time
    
    % Surface scaling constant
    sc = 5;
    
    % Surface curve
    f = x^2 + y^2;
    %f = cos(x);
    %f = 1/(x + y);
    
% -----------------------------------------------------------------------

% Find partial derivatives for gradient
dfdx = diff(f, x);
dfdy = diff(f, y);

% Create gradient function handles
dfdx_fun = matlabFunction(dfdx, 'Vars', [x, u2n, y, u4n]);
dfdy_fun = matlabFunction(dfdy, 'Vars', [x, u2n, y, u4n]);

% Create RK function handles
f1 = @(t, u1, u2, u3, u4) u2;
f2 = @(t, u1, u2, u3, u4) -g/m*dfdx_fun(u1, u2, u3, u4) - a/m*u2;
f3 = @(t, u1, u2, u3, u4) u4;
f4 = @(t, u1, u2, u3, u4) -g/m*dfdy_fun(u1, u2, u3, u4) - a/m*u4;

% Create surface function handle
f_fun = matlabFunction(f, 'Vars', [x, y]);

% Solve function
[t, x, ~, y, ~] = RK4_4(f1, f2, f3, f4, 0, xi, dxi, yi, dyi, h, t_f);

% Set up vectors and scalars
N = length(t);
z = zeros(1, N);

% Find trajectory
for i = 1:1:N
    z(i) = f_fun(x(i), y(i));
end

% Find surface
[xs, ys] = meshgrid(x(1:sc:end), y(1:sc:end));
zs = f_fun(xs, ys);

% Plot trajectory on surface
surf(xs, ys, zs);
colormap summer
shading interp
hold on;
plot3(x, y, z, 'r', 'LineWidth', 3);
hold off;

% Create motion plot
figure;
plot(x, y, 'LineWidth', 1);
grid on;
box on;

% Create trajectory plot
figure;
plot3(t, x, y, 'LineWidth', 1);
grid on;
box on;

% Create time plots
figure;
subplot(2, 1, 1);
plot(t, x, 'b', 'LineWidth', 1);
subplot(2, 1, 2);
plot(t, y, 'r', 'LineWidth', 1)

