function [K,F] = assemble1D(k, f, K, F, node_list)
    s = node_list(1);
    t = node_list(2);
    K(s:t, s:t) = K(s:t, s:t) + k;
    F(s:t) = F(s:t) + f;
end