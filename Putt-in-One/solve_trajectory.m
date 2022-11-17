% *********************************************************************** %
%                                                                         %
% Project           : Golf Trajectory Simulator                           %
%                                                                         %
% File name         : solve_trajectory.m                                  %
%                                                                         %
% Version           : 1.0                                                 %
%                                                                         %
% Author            : Gabriele J. Giuli                                   %
%                                                                         %
% Date created      : 28/07/2021                                          %
%                                                                         %
% *********************************************************************** %

function [t_de, x_de, y_de] = solve_trajectory(x, y, dfdx, dfdy, xi, dxi, yi, dyi, h, t_f, g, m, a)
    
    % Create RK function handles with linear interpolation
    f1 = @(t, u1, u2, u3, u4) u2;
    f2 = @(t, u1, u2, u3, u4) -g/m*interp2(x, y, dfdx, u1, u3) - a/m*u2;
    f3 = @(t, u1, u2, u3, u4) u4;
    f4 = @(t, u1, u2, u3, u4) -g/m*interp2(x, y, dfdy, u1, u3) - a/m*u4;
    
    % Solve differential equations
    [t_de, x_de, ~, y_de, ~] = RK4_4(f1, f2, f3, f4, 0, xi, dxi, yi, dyi, h, t_f);
    
end