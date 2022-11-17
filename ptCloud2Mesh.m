function [x_int, y_int, z_int] = ptCloud2Mesh(x_sct, y_sct, z_sct, x_step, y_step)
    
    % Create vectors
    x_int = min(x_sct):x_step:max(x_sct);
    y_int = min(y_sct):y_step:max(y_sct);
    
    % Find dimensions
    N = length(x_int);
    M = length(y_int);
    
    % Initialize output
    z_int = zeros(M, N);
    
    % Generate grid
    for i = 1:1:N
        for j = 1:1:M 
            z_int(j, i) = griddata(x_sct, y_sct, z_sct, x_int(i), y_int(j));
        end
    end
    
end