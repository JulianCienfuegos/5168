function [] = fem_1D_plot(X,U,psi)

%-----------------------------------------------------------------------
%
%  fem_1D_plot(X,U,psi)
%   
%  This function plots a 1D finite element solution. 
%
%  Written by Ethan Kubatko
%
%-----------------------------------------------------------------------
%
%  Input:
%  ------
%    X:   x coordinates of nodes
%    U:   finite element nodal solutions
%    psi: shape function structure as defined by lagrange_poly
%
%
%  Output:
%  -------
%    Produces the finite element
%
%-----------------------------------------------------------------------

% Determine degree of element
p = length(psi)-1;
% Create some points for plotting the solution
Xend = X(1:p:end);
nelems = (length(X)-1)/p;

% Divide the master element into enough points to get a smooth plot
m  = 20; xi = linspace(-1,1,m);

% Create figure
box on
hold on
ii = 1;
for j = 1:nelems
    Uh = zeros(1,p+1);
    for i = 1:p+1
        Uh = Uh + U(ii)*psi(i).fun;                
        ii = ii + 1;
    end    
    Uplot = polyval(Uh,xi);
    Xplot = linspace(Xend(j),Xend(j+1),m); 
    plot(Xplot,Uplot,'b','LineWidth',2)
    ii = ii-1;
end
plot(X,U,'bo','MarkerFaceColor','white','MarkerSize',10)
plot(Xend,U(1:p:end),'bo','MarkerFaceColor',[.49 1 .63],'MarkerSize',8)
xlabel('x','FontSize',13,'FontWeight','demi');
ylabel('Uh','FontSize',13,'FontWeight','demi')
axis([X(1), X(end), min(min(U)),1.2*max(max(U))]);
hold off
end
