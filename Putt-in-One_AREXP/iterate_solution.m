% *********************************************************************** %
%                                                                         %
% Project           : Golf Trajectory Simulator                           %
%                                                                         %
% File name         : iterate_solution.m                                  %
%                                                                         %
% Version           : 1.0                                                 %
%                                                                         %
% Author            : Gabriele J. Giuli                                   %
%                                                                         %
% Date created      : 28/07/2021                                          %
%                                                                         %
% *********************************************************************** %


function [throw_v_new, d, x_de, y_de] = iterate_solution(throw_v, cd_h, gamma, x_h, y_h, x, y, dfdx, dfdy, xi, yi, h, t_f, g, m, a)

    % Find distance at central point
    [~, x_c, y_c] = solve_trajectory(x, y, dfdx, dfdy, ...
            xi, throw_v(1), yi, throw_v(2), h, t_f, g, m, a);
    d_c = calculate_distance(x_c, y_c, x_h, y_h);
    
    % Find distance at x perturbed point
    [~, x_fx, y_fx] = solve_trajectory(x, y, dfdx, dfdy, xi, ...
        throw_v(1) + cd_h, yi, throw_v(2), h, t_f, g, m, a);
    d_fx = calculate_distance(x_fx, y_fx, x_h, y_h);
    
    % Find distance at y perturbed point
    [~, x_fy, y_fy] = solve_trajectory(x, y, dfdx, dfdy, xi, ...
        throw_v(1), yi, throw_v(2) + cd_h, h, t_f, g, m, a);
    d_fy = calculate_distance(x_fy, y_fy, x_h, y_h);
    
    % Find gradient components
    grad_local = [(d_fx - d_c) / cd_h; (d_fy - d_c) / cd_h];
    
    % Adjust guess solution according to gradient
    throw_v_new = throw_v - gamma*grad_local;
    
    % Assign previous output parameters
    d = d_c;
    x_de = x_c;
    y_de = y_c;
    
end