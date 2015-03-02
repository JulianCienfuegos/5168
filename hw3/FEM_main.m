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
    k_c = (2/h)*KofX(n)*ke.k + (h/2)*BofX(n)*ke.b;
    f_c = (h/2)*FofX(n)*fe;
    % Now plug these contributions into the appropriate place.
    [K, F] = assemble1D(k_c, f_c, K, F, node_list);
end
[K, F] = enforce_boundaries(K, F, NODEBC1, NODEBC2, VBC1, VBC2);
u = K\F;


d = length(u)
u_exact = @(x) log(1+x)/log(2) - x;
x = linspace(0, 1, d);
plot(x, u_exact(x), 'b')
hold on
plot(linspace(0, 1, d), u);

length(u)
length(u_exact(x))
disp(norm(u' - u_exact(x), inf))