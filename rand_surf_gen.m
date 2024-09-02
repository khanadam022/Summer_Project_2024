% given fixed gaussian acf and length scale, output/save a rand surf as symbolic fn
% this fn is called by averages.m
% We used this fn to create surf2.txt which is called by Numerical_method.m

function surf = rand_surf_gen()
    % increasing N really increases the jaggedness of the surface
    % don't increase length scale l beyond 5
    N = 25; delv = 0.01; l = 3; 
    v = [1:N]*delv; %nodes 
    
    % calculating the amplitudes 
    rho = @(x) exp(-x.^2/l^2).*cos(x*v);
    B = integral(rho, -inf, inf, 'ArrayValued',true);
    B = B*(2/pi);
    A = B.^0.5;
    max_amp = sum(A); 
    % random phases 
    Phase = rand(1, N)*2*pi; 
    
    % defining and adding the different sinusoidal components
    syms x; 
    h = cell(1,N);
    surf = 0;
    for n = 1:N
        h{n} = A(n)*sin(v(n)*x + Phase(n));
        surf = surf + h{n};
    end
    
    % get our surface bounded between +/- 0.5
    surf = surf * 0.5 / max_amp;
    
    % save to file
    % filename = 'surf2.txt'; 
    % save(filename, 'surf','-mat');
end 