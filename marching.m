% This program implements the MARCHING METHOD to reconstruct the surface at
% each node from the values scattered_field[1:500]
tic 
% define initial params
shift = 0.6; 
k = 1; w = 8; Z = 0.3+shift; Z0 = 20+shift; alpha =  0.5 * (1i/(2 * pi *k))^0.5; 
psi_inc = @(x,z) (w * (w^2 + 2i*x/k)^(-0.5) * exp( - (z - Z0)^2 / (w^2 + 2i*x/k))); 
N = 516; % for surf1.txt

% retrieve scattered field values
[~, scattered_field] = numerical_method(); 

% n = 1 values set to 0
normal_deriv = zeros(1,N); surf = zeros(1,N); 

N = 500; % max number of iterates 
for n = 2:N
    gamma = Gamma3([0, real(surf(1:n-1))], alpha, k, Z); 
    gamma = - gamma; % SIGN INVERSION
    sum = 0; 
    for r = 1:n-1
        sum = sum + gamma(r)*normal_deriv(r);
    end
    normal_deriv(n) = (scattered_field(n) - sum)/gamma(n);
    surf(n) = real(L_operator(normal_deriv(1:n), w, k, alpha, Z0, psi_inc)); % forcing this to be real
end 

% surf now contains calculated values
% calculate normal_deriv and surf again using these vals
% gamma3.m is called rather than Gamma3.m which sets h=0
% for n = 1:N
%     gamma = gamma3(real(surf(1:n)), alpha, k, Z); 
%     gamma = - gamma; % SIGN INVERSION
%     sum = 0; 
%     for r = 1:n-1
%         sum = sum + gamma(r)*normal_deriv(r);
%     end
%     normal_deriv(n) = (scattered_field(n) - sum)/gamma(n);
%     surf(n) = 1* real(L_operator(normal_deriv(1:n), w, k, alpha, Z0, psi_inc)); % scaled by 2
% end 

clf; 
plot(50:N, surf(50:N));
hold on; 
title('reconstruction for flat surface - 2 iterations');
%disp(real(surf(1:N)));

 % USING SURF1.TXT
fileID = fopen('surf1.txt');
surf = fscanf(fileID, '%f');
surf = surf *0.10; % scaling the surface
% PLOTTING THE SURFACE
plot(50:N, surf(50:N));
legend('Reconstruction', 'Actual');
toc