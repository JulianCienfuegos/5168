function [K, F] = BCTYPE1(nbc1, NODEBC1, VBC1, K, F, nnodes)
    for i = [1:nbc1]
        r = NODEBC1(i);
        K(r,:) = zeros(1,nnodes);
    end

    for k = 1:nbc1 
        i = NODEBC1(k);
        F = F - VBC1(k)*K(:,i);
    end

    for i = [1:nbc1]
        c = NODEBC1(i);
        K(:,c) = zeros(nnodes,1);
        K(c,c) = 1;
        F(c) = VBC1(i);
    end
end