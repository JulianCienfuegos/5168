function [ke, fe] = element1d(psi)
%{
	The shape functions will have degree d = length(psi) - 1.
	Therefore the product of two of these functions will have degree
	2*d.
	The derivatives of the shape functions will have degree d - 1. 
	Therefore the product of any two of these functions will have
	degree 2*d - 2. 
%}

d = length(psi) - 1;

n_kk = ceil((((2*d-2)+1) + 1)/2);
[x_gauss_kk, w_gauss_kk] = gauss_points_and_weights(n_kk);

n_kb = ceil((2*d+1)/2);
[x_gauss_kb, w_gauss_kb] = gauss_points_and_weights(n_kb);
 
n_f = ceil((d+1)/2);
[x_gauss_f, w_gauss_f] = gauss_points_and_weights(n_f);
ke.k = zeros(d+1);
ke.b = zeros(d+1);
fe = zeros(d+1, 1);

for i = 1:d+1
	for j = 1:i
		g = conv(psi(i).der, psi(j).der);
		g_gauss = polyval(g, x_gauss_kk);
		ke.k(i, j) = dot(w_gauss_kk, g_gauss);
		
		g = conv(psi(i).fun, psi(j).fun);
		g_gauss = polyval(g, x_gauss_kb);
		ke.b(i, j) = dot(w_gauss_kb, g_gauss);
	end
	for j = 1:i-1
		% Don't put this in the for loop above to avoid diagonal elements.
		ke.k(j, i) = ke.k(i, j);
		ke.b(j, i) = ke.b(i, j);
	end
	
	f_gauss = polyval(psi(i).fun, x_gauss_f);
	fe(i) = dot(w_gauss_f, f_gauss);
end

end
