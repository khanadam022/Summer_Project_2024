% Variant of Marching.m that calculates the average function distance
% between the actual and simulated surfaces across MAX = 100 iterations
tic
% define initial params
shift = 0.6;
k = 1; w = 8; Z = 0.5+shift; Z0 = 20+shift; alpha =  0.5 * (1i/(2 * pi *k))^0.5;
psi_inc = @(x,z) (w * (w^2 + 2i*x/k)^(-0.5) * exp( - (z - Z0)^2 / (w^2 + 2i*x/k)));

% this is what's changed
delx = 1.0;
N = floor(500/delx); 

total_distance = 0; 
total_acf_dist = 0; 
MAX = 2; % changed to a small number for testing purposes
for index = 1:MAX
    % retrieve scattered field values
    [surf2, ~, scattered_field] = Numerical_method_rand(delx);
    
    % n = 1 values set to 0 and perform marching alg
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
       %surf(n) = surf(n) + 0.1*(1/delx); % NO IDEA WHY THIS SHIFT WORKS BUT IT DOES
    end

    % plotting the surface to check
    %Nmin = floor(50/delx); 
    Nmin = 1; % over-writes previous value
   
    if index == 1
        clf;
        plot(delx*(Nmin:N), surf(Nmin:N)); %Nmin REPLACES 1
        hold on; 
        fplot(surf2, [0,500]); % 50 REPLACES 0
        legend('reconstructed', 'actual');
        title(['Surface reconstruction given delx=1 and Z=', num2str(Z-shift)]);
    end
    
    % calculating the L2 norm
    % sum = 0; 
    % for n = Nmin:N % Nmin replaces 1
    %     sum = sum + (surf(n) - surf2(n*delx))^2 * delx;
    % end
    % distance = sqrt(sum);
    % total_distance = total_distance + distance; 

    % calc normalized ACF for this iteration
    surf2_discretized = zeros(1, N);
    for n = 1:N
        surf2_discretized(n) = surf2(n*delx);
    end
    r = xcorr(surf, surf2_discretized, N, 'normalized');
    acf_dist = max(r);
    total_acf_dist = total_acf_dist + acf_dist;

    disp(['ITERATION ',num2str(index),' COMPLETED']);
end 
% compute average function distance for given value of delx
av_distance = total_distance/MAX;  
av_acf_dist = total_acf_dist/MAX;
%disp(['The average L2 dist given delx=', num2str(delx), ' over ', num2str(MAX), ' iterations is ', num2str(av_distance)]);
%disp(['The L2 dist over ', num2str(MAX),' iterations given delx=1 and Z=', num2str(Z-shift), ' is ', num2str(av_distance)]);
disp(['The ACF dist over ', num2str(MAX),' iterations given delx=1 and Z=', num2str(Z-shift), ' is ', num2str(av_acf_dist)]);
%disp(['The average ACF dist given delx=', num2str(delx), ' over ', num2str(MAX), ' iterations is ', num2str(av_acf_dist)]);
toc
