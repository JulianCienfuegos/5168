function [K, F] = BCTYPE1(nbc1, NODEBC1, VBC1, K, F, nnodes)
for i = NODEBC1 
    K(i,:) = zeros(nnodes);
    F = F - VBC1(i)*K(:,i);
end
end

function [K, F] = BCType2()
    % Do nothing!
end