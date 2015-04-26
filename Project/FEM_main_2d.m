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

for num = ['3']
    f = strcat('2D_Problems/Problem_2.1/Dirichlet/Mesh', num ,'/*');
    disp(f)
    copyfile(f, '.')

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
    U = exact_2D(XNODES,YNODES,'Problem 2.1 DBC ');
    figure()
    trimesh(CONN, XNODES, YNODES, U)
    plt2 = trimesh(CONN, XNODES, YNODES, u);
    ttl = strcat('Problem 2.1 Exact Mesh', num);
    title(ttl);
    name = strcat(ttl, '.jpg');
    saveas(plt2, name);

end