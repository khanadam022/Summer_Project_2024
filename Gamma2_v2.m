% VARIANT OF Gamma2.m called by Numerical_method.m
% key difference is the stepsize parameter delx

% returns a num_nodes x num_nodes matrix of coeffs
% used in the direct problem to calc the scattered field at height z=Z
% Gamma2(n,r) = int_{x_{r-1}}^{x_r} G[x_n, Z, x', h(x')]dx'
function gamma_matrix2 = Gamma2_v2(k, Z, alpha, num_nodes, delx)
    %define basic parameters 
    beta = 0.5 * k * Z^2 * ones(num_nodes,1); %516x1 - ALTERED
    beta_root = beta.^0.5; 
    beta_matrix = beta .* ones(1, num_nodes); % 516 identical col vectors
    beta_root_matrix = beta_root .* ones(1, num_nodes); % 516 col vectors
    b1 = zeros(num_nodes);
    b2 = zeros(num_nodes);
    for n = 1:num_nodes
        for r = 1:n
            b1(n,r) = ((n-r+1)*delx)^(-0.5); % delx ADDED
            if r<n
                b2(n,r) = ((n-r)*delx)^(-0.5); % delx ADDED
            end
        end
    end
    
    % calculate the fresnel part of the expression
    const = (2/pi)^0.5; 
    fres_arg1 = const *beta_root_matrix .* b1; 
    fres_arg2 = const *beta_root_matrix .* b2; 
    C1 = fresnelC(fres_arg1);
    C2 = fresnelC(fres_arg2);
    S1 = fresnelS(fres_arg1);
    S2 = fresnelS(fres_arg2);
    fres_expression = -4i * alpha * const^(-1) * beta_root_matrix .* (C2 + 1i*S2 - C1 - 1i*S1); 
    
    % calculate the exponential part of the expression 
    exp_arg1 = 1i*beta_matrix .* b1.^2; 
    exp_arg2 = 1i*beta_matrix .* b2.^2; 
    exp_expression = -2*alpha* (exp(exp_arg1)./b1 - exp(exp_arg2)./b2); 
    % combine all calculations 
    gamma_matrix2 = (exp_expression + fres_expression); 
    
    % correct the main diagonal - note: b1 consists of just 1's so ommit
    % make sure there are 0's above the main diagonal rather than NaN
    arg = const * beta_root; 
    fres_exp = -4i*alpha*const^(-1) * beta_root.*(0.5*(1+1i) - fresnelC(arg) - 1i*fresnelS(arg));
    exp_exp =  -2*alpha *exp(1i* beta); 
    diagonal = (exp_exp + fres_exp); 
    for n = 1:num_nodes
        gamma_matrix2(n,n) = diagonal(n); 
        for r = n:num_nodes
            if r>n
                gamma_matrix2(n,r) = 0; 
            end 
        end 
    end
    gamma_matrix2 = - gamma_matrix2; %SIGN CHANGE
end 