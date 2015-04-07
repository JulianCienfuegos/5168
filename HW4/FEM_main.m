close all, clear all;

% Set up mesh and input
generate_1D_mesh();
generate_1D_input();
read_1D_mesh();
read_1D_input();

K = zeros(nnodes);
F = zeros(nnodes, 1);
psi = lagrange_poly(p); 
[ke, fe] = element1d(psi);
for n = 1:nelems
    node_list = CONN(n,:);
    h = XNODES(node_list(2)) - XNODES(node_list(1)); 
    k = (2/h)*KofX(n)*ke.k + (h/2)*BofX(n)*ke.b;
    f = (h/2)*FofX(n)*fe;
    [K, F] = assemble1D(k, f, K, F, node_list);
end
[K, F] = enforce_boundaries(K, F, NODEBC1, NODEBC2, VBC1, VBC2, KofX);

% Calculate a solution
u = K\F;

% Exact Solution
w_0 = 38.4;
L = 10;
EI =  1;
deflection = 5*w_0*L^4/(384*EI);

% Compare the approximation to the analytical deflection.
approx = u(ceil(end/2));
disp('Percent error is :');
disp(100*abs(deflection - approx)/deflection);