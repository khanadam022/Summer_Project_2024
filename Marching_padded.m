% testing Numerical_method_rand_padded.m
tic 
% define initial params
shift = 0.6;
k = 1; w = 8; Z = 0.85+shift; Z0 = 20+shift; alpha =  0.5 * (1i/(2 * pi *k))^0.5;
psi_inc = @(x,z) (w * (w^2 + 2i*x/k)^(-0.5) * exp( - (z - Z0)^2 / (w^2 + 2i*x/k)));

% retrieving scattered field data - seems to be working properly
delx = 1.0;
[surf2, normal_deriv,scattered_field] = Numerical_method_rand_padded(delx);
N = max(size(normal_deriv));
Nmin = floor(50/delx);

% PERFORM THE MARCHING ALGORITHM
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
plot(delx*(1:N), surf2(1:N)-shift); % THIS SURF2 IS A VECTOR NOT A FN HANDLE
title('Padded surface reconstruction delx = 1.0, Z = 0.5')
toc

% calc normalized ACF for this iteration
N = floor(500/delx);
Nmin = floor(50/delx);
surf2_discretized = zeros(1, N);
for n = Nmin:N
    surf2_discretized(n) = surf2(n);
end
r = xcorr(surf(Nmin:N), surf2_discretized(Nmin:N), N, 'normalized');
acf_dist = max(r);
disp(acf_dist);