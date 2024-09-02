% THIS FUNCTION SOLVES THE DIRECT PROBLEM - VARIANT OF numerical_method.M
% takes a given step size delx
% uses surf2.txt sampled at delx*[1:N] where N = floor(500/delx)
%outputs normal deriv on the surface and scattered field at z=Z (sampled at the aforementionned nodes)

% we get instability for delx>5
function [normal_deriv,scattered_field] = Numerical_method(delx)
    % load surface and sample at delx*[1:N]
    file = which('surf2.txt');
    v = load(file, "-mat", "surf"); 
    surf2 = v.surf; % symbolic surface
    surf2 = 2*surf2; % rescaling
    surf2 = matlabFunction(surf2);

    num_nodes = floor(500/delx);
    surf = zeros(1, num_nodes);
    for x = 1:num_nodes
        surf(x) = surf2(x*delx);
    end 

    %z0: source height, Z: scattered data height, Z0: shift in Gaussian beam
    shift = 0.6; % ADDED TO AVOID NEGATIVE VALUES
    surf = surf +shift; 
    k = 1; w = 8; Z = 0.3+shift; Z0 = 20+shift; alpha =  0.5 * sqrt(1i/(2 * pi *k));
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
end