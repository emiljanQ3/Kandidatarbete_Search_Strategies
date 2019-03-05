
function [pos, rot, col] = moveAllAgents(pos,targetPos, rot_a ,obstacle, L, threshold, r_c, mapSize)

    col = zeros(3,1);
    numAgents = size(pos,1);
    targets = zeros(numAgents,2);
    
    for i = 1:numAgents
       [targets(i,:), rot(i)] = moveAgent(pos(i,:),targetPos(i,:), rot_a(i), obstacle, L, threshold);
    end
    
    % Detekterar just nu kolisioner genom hinder 
    
    [rot, col] = collision(targets,pos,rot,rot_a,numAgents,r_c);

    pos = targets;
end