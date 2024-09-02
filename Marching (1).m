% THIS IS THE BREAD AND BUTTER Marching.m PROGRAM

% this is the old marching.m just that we call Numerical_method_rand
% VARIANT OF marching.m designed to reconstruct the surface at arbitrary
% nodes delx*[1:N]
tic
% define initial params
shift = 0.6;
k = 1; w = 8; Z = 0.5+shift; Z0 = 20+shift; alpha =  0.5 * (1i/(2 * pi *k))^0.5;
psi_inc = @(x,z) (w * (w^2 + 2i*x/k)^(-0.5) * exp( - (z - Z0)^2 / (w^2 + 2i*x/k)));

% retrieve scattered field values
delx = 1.0;
N = floor(500/delx); % GENERALIZATION from numerical_method.m
[surf2, ~, scattered_field] = Numerical_method_rand(delx);

% n = 1 values set to 0
normal_deriv = zeros(1,N); surf = zeros(1,N);
for n = 2:N
  gamma = Gamma3_v2([0, real(surf(1:n-1))], alpha, k, Z, delx);
  gamma = - gamma; % SIGN INVERSION
  sum = 0;
  for r = 1:n-1
      sum = sum + gamma(r)*normal_deriv(r);
  end
  normal_deriv(n) = (scattered_field(n) - sum)/gamma(n);
  surf(n) = real(L_operator_v2(normal_deriv(1:n), w, k, alpha, Z0, psi_inc, delx)); % forcing this to be real
end

% examining reconstructed surface
clf;
plot(delx*(1:N), surf(1:N));
hold on; 
fplot(surf2, [0,500]);

% title('reconstruction for flat surface: delx = 2.0');
%
% % USING SURF2.TXT
% file = which('surf2.txt');
% v = load(file, "-mat", "surf");
% surf2 = v.surf; % symbolic surface
% surf2 = 2*surf2; % rescaling
% fplot(surf2, [0,500]);
% legend('Reconstruction', 'Actual');
% toc

