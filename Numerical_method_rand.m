% THIS FUNCTION SOLVES THE DIRECT PROBLEM - VARIANT OF Numerical_method.M
% takes a given step size delx
% RANDOMLY GENERATES a surface sampled at delx*[1:N] where N = floor(500/delx)
%outputs random surf as a matlab fn, normal deriv on the surface and scattered field at z=Z (sampled at the aforementionned nodes)

% called in Marching.m and Test.m
% we get instability for delx>5
function [surf2, normal_deriv,scattered_field] = Numerical_method_rand(delx)
    % load surface and sample at delx*[1:N]
    surf2 = rand_surf_gen();
    surf2 = matlabFunction(surf2);

    num_nodes = floor(500/delx);
    surf = zeros(1, num_nodes);
    for x = 1:num_nodes
        surf(x) = surf2(x*delx);
    end 

    %z0: source height, Z: scattered data height, Z0: shift in Gaussian beam
    shift = 0.6; % ADDED TO AVOID NEGATIVE VALUES
    surf = surf +shift; 
    k = 1; w = 8; Z = 0.5+shift; Z0 = 20+shift; alpha =  0.5 * sqrt(1i/(2 * pi *k));
    psi_inc = @(x,z) (w * (w^2 + 2i*x/k)^(-0.5) * exp( - (z - Z0)^2 / (w^2 + 2i*x/k))); 
    
    % populate vector of Psi_inc - working + fast
    Psi_inc = zeros(num_nodes, 1); 
    for n = 1:num_nodes
        Psi_inc(n) = psi_inc(n*delx, surf(n)); % MODIFIED FROM numerical_method.m
    end
    
    % CALCULATES THE NORMAL DERIVATIVES
    gamma_matrix = Gamma_v2(k, surf, alpha, num_nodes, delx); 
    normal_deriv = - gamma_matrix\Psi_inc; % NEGATIVE SIGN ADDED

    %CALCULATE THE SCATTERED FIELD VALS AT z=Z
    gamma_matrix2 = Gamma2_v2(k, Z, alpha, num_nodes, delx);
    scattered_field = - gamma_matrix2 * normal_deriv; % NEGATIVE SIGN ADDED
    
    % plot(delx*(1:num_nodes), abs(scattered_field));
    % hold on; 
end