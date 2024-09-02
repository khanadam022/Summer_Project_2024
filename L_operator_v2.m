%VARIANT OF L_operator.m called in Marching.m
% now delx is an input parameter 

% given [normal_deriv(1), ..., normal_deriv(n)], calculates h_n = h(x_n)
% this fn is called in the marching.m problem
% this function takes col vectors
function [h_n, root_difference, sum] = L_operator_v2(normal_deriv, w, k, alpha, Z0, psi_inc, delx)
    sz = size(normal_deriv);
    n = max(sz);

    root_difference = zeros(1,n); 
    sum = 0;
    for r = 1: n
        root_difference(r) = sqrt((n-r+1)*delx) - sqrt((n-r)*delx); %MODIFIED to reflect sqrt(x_n - x_{r-1}) - sqrt(x_n - x_r)
        sum = sum + normal_deriv(r) .* root_difference(r); 
    end 

    factor = (w^2 + 2i*n*delx/k)/(2*Z0); % MODIFIED 
    h_n = (-2*alpha * sum/psi_inc(n*delx,0) - 1) * factor; % MODIFIED

    % mysterious rescaling + shift
    shift = 0.6; 
    h_n = -  0.5*(w^2/Z0 + (h_n)) - shift;
end 