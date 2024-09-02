% Takes a vector of height values h(x_1), ..., h(x_{n})
% Calcs the vector gamma_{n,0}, ..., gamma_{n,r}, ..., gamma_{n.n-1} used to
% calc Psi_prime(x_n)
% This program is called in marching.m to solve the inverse problem

function gamma = Gamma3(h, alpha, k, Z)
    sz = size(h);
    n = max(sz);
    h = zeros(1,n); %ADDED AS AN EXPERIMENT - SOMEHOW SOLVES EVERYTHING
    
    % gamma_{1,1} case 
    if n == 1
        %b1 = 1 ommitted from formula
        a = 0; 
        b = 0.5*k * Z^2; % assume h(0) = h_prime(0) = 0
        f1 =  2*alpha *exp( - 0.5*1i*k*a); % ALPHA AND 1i ADDED
        f2 = 2i*sqrt(b * pi/2);
        arg1 = sqrt(b * 2/pi);
        C1 = fresnelC(real(arg1));
        S1 = fresnelS(real(arg1));
       
        % short form of the expression 
        gamma = f1 * (exp(1i * b) + f2 * (1 * 0.5*(1+1i) - C1 - 1i*S1)); %sqrt(2/pi) replaced by 1
    else % n>1
        for j = 0:n-1
            % gamma_{n,1}
            if j == 0 
                % initial params
                b1 = 1/n; 
                b2 = 1/(n-1);
                a = 0;
                b = 0.5*k*Z^2; % assume h(0) = h_prime(0) = 0

                f1 = 2*alpha*exp( - 0.5*1i*k*a); %ALPHA AND 1i ADDED
                f2 = 2i*sqrt(b*pi/2);
                arg2 = sqrt(b * b2 * 2/pi); 
                arg1 = sqrt(b * b1 * 2/pi); 
                C2 = fresnelC(real(arg2)); 
                C1 = fresnelC(real(arg1));
                S2 = fresnelS(real(arg2)); 
                S1 = fresnelS(real(arg1));

                % long form of the expression
                gamma(1) = f1 * (exp(1i* b *b1)/sqrt(b1) - exp(1i* b *b2)/sqrt(b2)  + f2 * (C2 + 1i*S2 - C1 - 1i*S1));

            %gamma_{n,n}
            elseif j == n-1 
                % define the initial params
                h_prime = h(n) - h(n-1);
                %b1 = 1 ommitted from the formula
                a = -2 * h_prime * (Z - h(n));
                b = 0.5*k* (Z - h(n)) * (Z - h(n) - 2*h_prime * (1));

                f1 =  2*alpha *exp( - 0.5*1i*k*a);
                f2 = 2i*sqrt(b * pi/2);
                arg1 = sqrt(b * 2/pi);
                C1 = fresnelC(real(arg1));
                S1 = fresnelS(real(arg1));

                % short form of the expression 
                gamma(n) = f1 * (exp(1i*b) + f2 * (1 * 0.5*(1+1i) - C1 - 1i*S1)); % sqrt(2/pi) replaced by 1

            % THIS IS WHERE THE BULK OF THE CALCS ARE PERFORMED
            %calcs gamma_{n,j} where 0<j <n-1 
            else 
                % define initial params 
                h_prime = h(j+1) - h(j);
                b1 = 1/ (n-j); 
                b2 = 1/(n-j-1);
                a = -2 * h_prime * (Z - h(j));
                b = 0.5*k* (Z - h(j)) * (Z - h(j) - 2* h_prime * (n-j));
                
                f1 = 2*alpha*exp( - 0.5*1i*k*a);
                f2 = 2i*sqrt(b*pi/2);
                arg2 = sqrt(b * b2 * 2/pi); 
                arg1 = sqrt(b * b1 * 2/pi); 
                C2 = fresnelC(real(arg2)); 
                C1 = fresnelC(real(arg1));
                S2 = fresnelS(real(arg2)); 
                S1 = fresnelS(real(arg1));

                % long version of the formula
                gamma(j+1) = f1 * (exp(1i * b * b1)/sqrt(b1) - exp(1i * b * b2)/sqrt(b2)  + f2 * (C2 + 1i*S2 - C1 - 1i*S1));
            end
        end
    end
end