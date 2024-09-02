% THIS FUNCTION SOLVES THE DIRECT PROBLEM
% uses surf1.txt sampled at 1:516 
% outputs normal deriv on the surface and scattered field at z=Z (also sampled @ 1:516)
% SHIFTS ADDED

% called in marching.m
function [normal_deriv, scattered_field] = numerical_method()
    % USING SURF1.TXT
    fileID = fopen('surf1.txt');
    surf = fscanf(fileID, '%f');
    surf = surf * 0.1; % scaling the surface
    sz = size(surf);
    num_nodes = max(sz);
    % PLOTTING THE SURFACE
    % xx = 1:15;
    % plot(xx, surf(1:15));
    
    %z0: source height, Z: scattered data height, Z0: shift in Gaussian beam
    shift = 0.6; % ADDED TO AVOID NEGATIVE VALUES
    surf = surf +shift; 
    k = 1; w = 8; Z = 0.3+shift; Z0 = 20+shift; alpha =  0.5 * sqrt(1i/(2 * pi *k));
    psi_inc = @(x,z) (w * (w^2 + 2i*x/k)^(-0.5) * exp( - (z - Z0)^2 / (w^2 + 2i*x/k))); 
        
    % populate vector of Psi_inc - working + fast
    Psi_inc = zeros(num_nodes, 1); 
    for n = 1:num_nodes
        Psi_inc(n) = psi_inc(n, surf(n));
    end
    
    % CALCULATES THE NORMAL DERIVATIVES
    gamma_matrix = Gamma(k, surf, alpha, num_nodes); 
    normal_deriv = - gamma_matrix\Psi_inc; % NEGATIVE SIGN ADDED
    
    %CALCULATE THE SCATTERED FIELD VALS AT z=Z
    gamma_matrix2 = Gamma2(k, Z, alpha, num_nodes);
    scattered_field = - gamma_matrix2 * normal_deriv; % NEGATIVE SIGN ADDED

    % write the scattered field data to file 
    %writematrix(scattered_field, 'scattered_field_data.txt');
    %writematrix(normal_deriv, 'normal_deriv_data.txt');
end 

% TESTING REFLECTION PROPERTY
% xx = 1:num_nodes; 
% psi_inc_reflected = zeros(1, num_nodes);
% for n = 1:num_nodes
%     psi_inc_reflected(n) = - psi_inc(n, -Z);
% end 
% 
% clf; 
% plot(xx, real(scattered_field(1: num_nodes))); 
% hold on; 
% plot(xx, real(psi_inc_reflected(1:num_nodes)));
% title('Re pt of scattered field and psi inc reflected with Z = 0.3');


% TESTING NORMAL DERIVATIVE ON SURFACE PROPERTY
% clf; 
% xx = 1:num_nodes; 
% psi_inc_deriv_surface = zeros(1, num_nodes);
% for n = 1:num_nodes
%     psi_inc_deriv_surface(n) = psi_inc_deriv(n,0);
% end
% plot(xx, real(normal_deriv));
% hold on; 
% plot(xx, real(2*psi_inc_deriv_surface));
% title('Surface derivative Condition - Re part');