function [K, F] = enforce_boundaries(K, F, NODEBC1, NODEBC2, VBC1, VBC2, KofX)
    for i = 1:length(NODEBC1)
        idx = NODEBC1(i);             % Where is the dirichlet condition?
        if(idx~=0)		
            F = F - VBC1(i)*K(:,idx); % Modify F
            F(idx) = VBC1(i);
            K(idx,:) = 0;             % Zero out the row and column.
            K(:,idx) = 0;
            K(idx, idx) = 1;          % Set diagonal entry to one.
        end
    end
    for i = 1:length(NODEBC2)
        idx = NODEBC2(i);             % where is the neumann condition?
        if (idx == 1)                 % modify f
            F(1) = F(1) - KofX(1)*VBC2(i);
        elseif (idx > 1)
            F(end) = F(end) + KofX(end)*VBC2(i);
        end        
    end
end
