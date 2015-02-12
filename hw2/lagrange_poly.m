function psi = lagrange_poly(k)
%{
    k is the degree of the basis functions.
    First, set up the xi and y coordinates for the polyfit.
    Then, polyfit the coordinates for each set of y coordinates
    and this is a psi function. Use polyder to get the 
    derivtives.
%}

xi = linspace(-1, 1, k+1); 
for func=1:k+1
    y = zeros(1,k+1);
    y(func) = 1;
    F = polyfit(xi, y, k);
    psi(func).fun = F;
    psi(func).der = polyder(F);

end