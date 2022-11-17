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
    
    % Solving parameters
    h = 1E-2;   % Step size
    tau = 2;    % Number of time contants
    
% -----------------------------------------------------------------------

% Generate mesh
x = x_i:x_s:x_f;
y = y_i:y_s:y_f;

% Find values of surface
[x_mtx, y_mtx] = meshgrid(x, y);
z_mtx = f(x_mtx, y_mtx);

% Find numerical gradient and final execution time
[dfdx, dfdy] = gradient(z_mtx, x, y);
t_f = tau*m/a;

[a, b, c] = solve_trajectory(x, y, dfdx, dfdy, 0.9, -20, 0.9, -20, h, t_f, g, m, a);