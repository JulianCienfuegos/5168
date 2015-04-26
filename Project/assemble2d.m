function [K,F] = assemble2d(nodelist, ke, fe, K, F)
for r = [1:3]
    i = nodelist(r);
    F(i) = F(i) + fe(r);
    for s = [1:3]
        j = nodelist(s);
        K(i,j) = K(i,j) + ke(r,s);
    end
end
end