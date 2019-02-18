function pos = moveAllAgents(pos,targetPos, obstacle, L, threshold, r_c)
    
    numAgents = size(pos,1);
    targets = zeros(numAgents,2);
    
    for i = 1:numAgents
       targets(i,:) = moveAgent(pos(i,:),targetPos(i,:), obstacle,L, threshold);
    end
    
    for i = 1:numAgents-1
        for j = i+1:numAgents
            if norm(targets(i) - targets(j)) < r_c
                pos_1 = targets(i,:);
                scatter(pos_1(1), pos_1(2), 'r');
            end
        end
    end
    
    pos = targets;
end