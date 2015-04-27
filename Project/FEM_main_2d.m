%
% This is my 2d FEM solver.
%
% By Melvyn Ian Drag
% 26 April 2015
%

close all; clear all;
set(0,'DefaultAxesFontSize', 18)

L2 = [];
Linf = []; 
% % 
for num = [ '1', '2', '3']
    f = strcat('2D_Problems/Problem_2.1/Dirichlet/Mesh', num ,'/*');
    disp(f)
    copyfile(f, '.')
    probid = 'Problem 2.1 DBC ';    
    
    read_2D_mesh();
    read_2D_input();
    K = zeros(nnodes);
    F = zeros(nnodes,1);

    for n=1:nelems
        nodelist = CONN(n,:);
        x = XNODES(nodelist);
        y = YNODES(nodelist);
        k = KofXY(n);
        b = BofXY(n);
        f = FofXY(n);
        [ke, fe] = element2d(x, y, k, b, f);
        [K,F] = assemble2d(nodelist, ke, fe, K, F);
    end

    [K, F] = BCTYPE1(nbc1, NODEBC1, VBC1, K, F, nnodes);
    [K, F] = BCType2(K, F);

    % Compute and plot numerical solution ---------------------------------
    u = K\F;
    figure()
    plt = trimesh(CONN, XNODES, YNODES, u);
    ttl = strcat('Problem 2.1 Mesh', num);
    title(ttl);
    name = strcat(ttl, '.jpg');
    saveas(plt, name)
    
    % Compute and plot exact solution -------------------------------------
    U = exact_2D(XNODES,YNODES,probid);
    figure()
    trimesh(CONN, XNODES, YNODES, U)
    plt2 = trimesh(CONN, XNODES, YNODES, u);
    ttl = strcat('Problem 2.1 Exact Mesh', num);
    title(ttl);
    name = strcat(ttl, '.jpg');
    saveas(plt2, name);
    
    % Compute error norms -------------------------------------------------
    L2 = 0;
    Linf = 0;
    
    for j = 1:nelems
        nodelist = CONN(j,:);
        xi = XNODES(nodelist);
        yi = YNODES(nodelist);
        ui = u(nodelist);
       
        A_e = 0.5*(xi(1)*yi(2) - yi(1)*xi(2) + xi(2)*yi(3) - yi(2)*xi(3) + ...
            xi(3)*yi(1) - yi(3)*xi(1));
        w = [1/6, 1/6, 1/6];
        p1 = @(ksi, eta) 1 - ksi - eta;
        p2 = @(ksi, eta) ksi;
        p3 = @(ksi, eta) eta;
        psi = [p1(0.5, 0), p1(0.5, 0.5), p1(0, 0.5);...
               p2(0.5, 0), p2(0.5, 0.5), p2(0, 0.5);...
               p3(0.5, 0), p3(0.5, 0.5), p3(0, 0.5)];
        for k = [1:3]
            u_h = dot(ui, psi(:,k));
            x = dot(xi, psi(:,k));
            y = dot(yi, psi(:,k));
            
            U = exact_2D(x, y, probid);
            
            L2 = L2 + w(k) * 2*A_e*(U - u_h)^2;
            Linf = max(Linf, abs(U - u_h));
        end
    end
    L2 = sqrt(L2);
end