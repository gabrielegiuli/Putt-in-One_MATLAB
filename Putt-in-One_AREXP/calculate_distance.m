% *********************************************************************** %
%                                                                         %
% Project           : Golf Trajectory Simulator                           %
%                                                                         %
% File name         : calculate_distance.m                                %
%                                                                         %
% Version           : 1.0                                                 %
%                                                                         %
% Author            : Gabriele J. Giuli                                   %
%                                                                         %
% Date created      : 28/07/2021                                          %
%                                                                         %
% *********************************************************************** %

function d = calculate_distance(x, y, x_h, y_h)

    % Find trajectory distances from hole
    ds = sqrt((x - x_h).^2 + (y - y_h).^2);
    
    % Find indexes of closest point
    [~, index_1] = min(ds);
    
    % Remove closest point from array and find second to closest point
    ds(index_1) = max(ds) + 1;
    [~, index_2] = min(ds);
    
    % Extract point coordinates
    x_1 = x(index_1);
    y_1 = y(index_1);
    x_2 = x(index_2);
    y_2 = y(index_2);
    
    % Find line passing through the closest points
    m = (y_2 - y_1)/(x_2 - x_1);
    q = y_1 - m*x_1;
    
    % Perform distance linear interpolation
    d = abs(y_h - (m*x_h + q))/sqrt(1 + m^2);
    
end