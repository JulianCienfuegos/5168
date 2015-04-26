close all, clear all;
set(0,'DefaultAxesFontSize', 18)

dir = ['1', '2', '3', '4', '5', '6'];
err_l = [];
err_r = [];
H = [];

for d = dir
    f = strcat('1D_Problems/Problem_1.2/Mesh_',d,'/*');
    disp(f)
    copyfile(f, '.')

    % Set up mesh and input
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

    ui = 100;
    uo = 0;
    ri = 1;
    ro = 10;
    sig_exact = 1.7*(ui-uo)/log(ro/ri);
    sig = @(r, u_prime) -1.7*r*u_prime;
    u_pl = (u(2) - u(1)) / h;
    u_pr = (u(end) - u(end-1)) / h;
    sig_l = sig( 1, u_pl);
    sig_r = sig(10, u_pr);

    rel_err_l = (abs((sig_exact - sig_l)/sig_exact));
    rel_err_r = (abs((sig_exact - sig_r)/sig_exact));

    err_l = [err_l, rel_err_l];
    err_r = [err_r, rel_err_r];
    H = [H, h];
end

log_err_l = fliplr(log(err_l));
log_err_r = fliplr(log(err_r));
log_H     = fliplr(log(H));
pl    = polyfit(log_H, log_err_l, 1);
pr    = polyfit(log_H, log_err_r, 1);
plv   = polyval(pl, log_H);
prv   = polyval(pr, log_H);

figure()
plot(log_H, log_err_l, log_H, plv)
title('Convergence at Left Node')
xlabel('log(h)')
ylabel('log(error)')
legend('Error', 'Polyfit')

figure()
plot(log_H, log_err_r, log_H, prv)
title('Convergence at Right Node')
xlabel('log(h)')
ylabel('log(error)')
legend('Error', 'Polyfit')

disp(pl(1))
disp(pr(1))