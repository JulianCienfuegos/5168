function psi = lagrange_poly2(k)
%{
k is the degree of the basis functions.
The general structure of this code is:
1) Set some necessary variables.
2) To make the coefficients we first convolute the numerator terms in the 
   expression for the lagrange polynomial.
   Then, we multiply the denominator terms and multiply the denominator
   by the coefficient array. 
   We have to go through this for each basis function.
3) Store the basis functions and their derivatives.
%}

xi_m = -1; % minimum xi
xi_M =  1; % maximum xi
xi = linspace(xi_m, xi_M, k+1); % generate k+1 nodes.

for func=1:k+1
    
    % Initialize identity elements.
    f = [1];
    denom = 1;
    
    % In this loop we will multiply the numerator and denominator factors.
    for factor = 1:k+1
        if (factor ~= func)
            f = conv(f, [1, -xi(factor)]);
            denom = denom * (xi(func) - xi(factor)); 
        end
    end
    
    % Set the basis function and its derivative.
    F = f/denom;
    psi(func).fun = F;
    power_rule = [k:-1:0];
    dF = F.*power_rule;          % Apply power rule.
    psi(func).der = dF(1:end-1); % Trim off end 0 element.

end