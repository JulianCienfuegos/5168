close all, clear all;
%{
    Solve a problem using the FEM.
%}
read_1D_mesh();
read_1D_input();
K = zeros(nnodes);
F = zeros(nnodes, 1);
psi = lagrange_poly(p); % p is assigned in read_1d_mesh
[ke, fe] = element1d(psi);
for n = 1:nelems
    node_list = CONN(n,:);
    h = XNODES(node_list(2)) - XNODES(node_list(1)); 
    % Contributions to the global matrix
    k = (2/h)*KofX(n)*ke.k + (h/2)*BofX(n)*ke.b;
    f = (h/2)*FofX(n)*fe;
    % Now plug these contributions into the appropriate place.
    [K, F] = assemble1D(k, f, K, F, node_list);
end
[K, F] = enforce_boundaries(K, F, NODEBC1, NODEBC2, VBC1, VBC2, KofX);

% Calculate a solution
u = K\F;

% Compare the solution to the real thing. View the plot and the error.
d = length(u);
u_exact = @(x) log(1+x)/log(2) - x;
%u_exact = @(x) x + 2 + sin(x)*sec(1);
x = linspace(0, 1, d);
figure()
plot(x, u_exact(x), 'b', x, u, 'r')
legend('y = u\_exact(x)','y = u\_approx','Location','northeast')
title('Comparison of FEM Solution to Exact Solution.')
x1 = 0.7;
y1 = 0.075;
str1 = strcat('Norm of error: ', num2str(norm(u' - u_exact(x), inf)));
text(x1, y1, str1);
xlabel('x')
ylabel('u(x)')

norm(u' - u_exact(x), inf)