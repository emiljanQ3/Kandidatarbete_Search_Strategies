function [pos rot_a col] = moveAllAgents(pos,targetPos, rot_a ,obstacle, L, threshold, r_c)
    col = zeros(3,1);
    numAgents = size(pos,1);
    targets = zeros(numAgents,2);
    
    for i = 1:numAgents
       targets(i,:) = moveAgent(pos(i,:),targetPos(i,:), obstacle,L, threshold);
    end
    
    new_rot_a = rot_a;
    k = 0;
    for i = 1:numAgents-1
        for j = i+1:numAgents
            if norm(targets(i,:) - targets(j,:)) < r_c  
                k = k+1;
                col(1:2,1) = (targets(i,:) + targets(j,:))/2;
                col(3,1) = 1;
                
                new_rot_a(j) = 2*pi*rand;
                new_rot_a(i) = 2*pi*rand;
                targets(j,:) = pos(j,:);
                targets(i,:) = pos(i,:);
            end
        end
    end
    rot_a = new_rot_a;
    pos = targets;
end