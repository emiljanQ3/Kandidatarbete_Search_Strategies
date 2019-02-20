function [pos, rot_a] = moveAllAgents(pos,targetPos, rot_a ,obstacle, L, threshold, r_c)
    
    numAgents = size(pos,1);
    targets = zeros(numAgents,2);
    
    for i = 1:numAgents
       [targets(i,:), rot_a(i)] = moveAgent(pos(i,:),targetPos(i,:), rot_a(i), obstacle, L, threshold);
    end
    
    
    for i = 1:numAgents-1
        for j = i+1:numAgents
            if norm(targets(i,:) - targets(j,:)) < r_c                
                rot_a(j) = 2*pi*rand;
                rot_a(i) = 2*pi*rand;
                targets(j,:) = pos(j,:);
                targets(i,:) = pos(i,:);                
            end
        end
    end
    pos = targets;
end