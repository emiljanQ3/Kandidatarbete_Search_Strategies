
function [pos, rot, col] = moveAllAgentsCircle(pos,targetPos, rot_a , R, r_c)

    numAgents = size(pos,1);
    targets = zeros(numAgents,2);
    
    for i = 1:numAgents
       [targets(i,:), rot(i)] = moveAgentCircle(pos(i,:),targetPos(i,:), rot_a(i), R);
    end
    
    % Detekterar just nu kolisioner genom hinder 
    
    [rot, col] = collision(targets,pos,rot,rot_a,numAgents,r_c);

    pos = targets;
end