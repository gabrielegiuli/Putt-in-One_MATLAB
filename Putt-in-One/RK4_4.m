% *********************************************************************** %
%                                                                         %
% Project           : Golf Trajectory Simulator                           %
%                                                                         %
% File name         : RK4_4.m                                             %
%                                                                         %
% Version           : 1.0                                                 %
%                                                                         %
% Author            : Gabriele J. Giuli                                   %
%                                                                         %
% Date created      : 19/07/2021                                          %
%                                                                         %
% *********************************************************************** %

function [t, u1, u2, u3, u4] = RK4_4(f1, f2, f3, f4, ti, u1i, u2i, u3i, u4i, h, tf)
    
    % Determine number of steps N
    N = round((tf - ti)/h);
    
    % Initialize vectors
    t = [ti, zeros(1, N-1)];
    u1 = [u1i, zeros(1, N-1)];
    u2 = [u2i, zeros(1, N-1)];
    u3 = [u3i, zeros(1, N-1)];
    u4 = [u4i, zeros(1, N-1)];
    
    % Apply RK to calculate x_(i+1), y_(i+1) and z_(i+1)
    for i = 2:N
        k1_u1 = f1(t(i-1), u1(i-1), u2(i-1), u3(i-1), u4(i-1));
        k1_u2 = f2(t(i-1), u1(i-1), u2(i-1), u3(i-1), u4(i-1));
        k1_u3 = f3(t(i-1), u1(i-1), u2(i-1), u3(i-1), u4(i-1));
        k1_u4 = f4(t(i-1), u1(i-1), u2(i-1), u3(i-1), u4(i-1));
        
        k2_u1 = f1(t(i-1) + 0.5*h, u1(i-1) + 0.5*k1_u1*h, u2(i-1) + 0.5*k1_u2*h, u3(i-1) + 0.5*k1_u3*h, u4(i-1) + 0.5*k1_u4*h);
        k2_u2 = f2(t(i-1) + 0.5*h, u1(i-1) + 0.5*k1_u1*h, u2(i-1) + 0.5*k1_u2*h, u3(i-1) + 0.5*k1_u3*h, u4(i-1) + 0.5*k1_u4*h);
        k2_u3 = f3(t(i-1) + 0.5*h, u1(i-1) + 0.5*k1_u1*h, u2(i-1) + 0.5*k1_u2*h, u3(i-1) + 0.5*k1_u3*h, u4(i-1) + 0.5*k1_u4*h);
        k2_u4 = f4(t(i-1) + 0.5*h, u1(i-1) + 0.5*k1_u1*h, u2(i-1) + 0.5*k1_u2*h, u3(i-1) + 0.5*k1_u3*h, u4(i-1) + 0.5*k1_u4*h);
        
        k3_u1 = f1(t(i-1) + 0.5*h, u1(i-1) + 0.5*k2_u1*h, u2(i-1) + 0.5*k2_u2*h, u3(i-1) + 0.5*k2_u3*h, u4(i-1) + 0.5*k2_u4*h);
        k3_u2 = f2(t(i-1) + 0.5*h, u1(i-1) + 0.5*k2_u1*h, u2(i-1) + 0.5*k2_u2*h, u3(i-1) + 0.5*k2_u3*h, u4(i-1) + 0.5*k2_u4*h);
        k3_u3 = f3(t(i-1) + 0.5*h, u1(i-1) + 0.5*k2_u1*h, u2(i-1) + 0.5*k2_u2*h, u3(i-1) + 0.5*k2_u3*h, u4(i-1) + 0.5*k2_u4*h);
        k3_u4 = f4(t(i-1) + 0.5*h, u1(i-1) + 0.5*k2_u1*h, u2(i-1) + 0.5*k2_u2*h, u3(i-1) + 0.5*k2_u3*h, u4(i-1) + 0.5*k2_u4*h);
        
        k4_u1 = f1(t(i-1) + h, u1(i-1) + k3_u1*h, u2(i-1) + k3_u2*h, u3(i-1) + k3_u3*h, u4(i-1) + k3_u4*h);
        k4_u2 = f2(t(i-1) + h, u1(i-1) + k3_u1*h, u2(i-1) + k3_u2*h, u3(i-1) + k3_u3*h, u4(i-1) + k3_u4*h);
        k4_u3 = f3(t(i-1) + h, u1(i-1) + k3_u1*h, u2(i-1) + k3_u2*h, u3(i-1) + k3_u3*h, u4(i-1) + k3_u4*h);
        k4_u4 = f4(t(i-1) + h, u1(i-1) + k3_u1*h, u2(i-1) + k3_u2*h, u3(i-1) + k3_u3*h, u4(i-1) + k3_u4*h);
        
        u1(i) = u1(i-1) + h*(1/6*k1_u1 + 1/3*k2_u1 + 1/3*k3_u1 + 1/6*k4_u1);
        u2(i) = u2(i-1) + h*(1/6*k1_u2 + 1/3*k2_u2 + 1/3*k3_u2 + 1/6*k4_u2);
        u3(i) = u3(i-1) + h*(1/6*k1_u3 + 1/3*k2_u3 + 1/3*k3_u3 + 1/6*k4_u3);
        u4(i) = u4(i-1) + h*(1/6*k1_u4 + 1/3*k2_u4 + 1/3*k3_u4 + 1/6*k4_u4);
        
        t(i) = t(i-1) + h;
        
        % Check for NaN
        if any(isnan([u1(i), u2(i), u3(i), u4(i)]))
            
            % Delete remaining part of vectors
            u1(i:end) = [];
            u2(i:end) = [];
            u3(i:end) = [];
            u4(i:end) = [];
            
            % Exit function
            return;
        end
        
    end
end