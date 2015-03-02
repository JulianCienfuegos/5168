function [K, F] = enforce_boundaries(K, F, NODEBC1, NODEBC2, VBC1, VBC2)
if ~isempty(NODEBC1)
    for i = 1:length(NODEBC1)
        % Where is the dirichlet condition?
        idx = NODEBC1(i);
        
                
        % Modify F
        F = F - VBC1(i)*K(:,idx);
        F(idx) = VBC1(i);
        

        % Zero out the row and column.
        K(idx,:) = 0;
        K(:,idx) = 0;
        
        % Set diagonal entry to one.
        K(idx, idx) = 1;
    end
end

if ~isempty(NODEBC2)
% Do other stuff
end

end