
function [pos, rot_a, col] = moveAllAgents(pos,targetPos, rot_a ,obstacle, L, threshold, r_c, mapSize)

    col = zeros(3,1);
    numAgents = size(pos,1);
    targets = zeros(numAgents,2);
    
    for i = 1:numAgents
       [targets(i,:), rot_a(i)] = moveAgent(pos(i,:),targetPos(i,:), rot_a(i), obstacle, L, threshold);
    end
    
    % Detekterar just nu kolisioner genom hinder 
    for i = 1:numAgents-1
        for j = i+1:numAgents
            if norm(targets(i,:) - targets(j,:)) < r_c    
                col(1:2) = (targets(i,:)+targets(j,:))/2;
                col(3) = 1;
                rot_a(j) = 2*pi*rand;
                rot_a(i) = 2*pi*rand;
            end
        end
    end


    pos = targets;
end