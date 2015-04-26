function [ke, fe] = element2d(x, y, k, b, f)
    A_e = 0.5*(x(1)*y(2) - y(1)*x(2) + x(2)*y(3) - y(2)*x(3) + ...
        x(3)*y(1) - y(3)*x(1));% Compute area
    
    X = constructMat(x);
    Y = constructMat(y);
    
    w = [1/6, 1/6, 1/6];
    p1 = @(ksi, eta) 1 - ksi - eta;
    p2 = @(ksi, eta) ksi;
    p3 = @(ksi, eta) eta;
    psi = [p1(0.5, 0), p1(0.5, 0.5), p1(0, 0.5);...
           p2(0.5, 0), p2(0.5, 0.5), p2(0, 0.5);...
           p3(0.5, 0), p3(0.5, 0.5), p3(0, 0.5)];
       
    for i = [1:3]
        for j = [1:3]
            g_gauss = psi(i,:).*psi(j,:);
            ke(i,j) = k*(X(i,j) + Y(i,j))/(4*A_e) + ...
                2*b*A_e*dot(g_gauss, w);
        end
        fe(i) = f*A_e/3;
    end
end

function M = constructMat(m)
    M = zeros(3);
    M(1,1) = (m(2) - m(3))^2;
    M(1,2) = (m(2) - m(3))*(m(3) - m(1));
    M(2,1) = M(1,2);
    M(1,3) = (m(2) - m(3))*(m(1) - m(2));
    M(3,1) = M(1,3);
    M(2,2) = (m(3) - m(1))^2;
    M(2,3) = (m(3) - m(1))*(m(1) - m(2));
    M(3,2) = M(2,3);
    M(3,3) = (m(1) - m(2))^2;
end