function [rot_a,col] = collision(targets,pos,rot,rot_a,numAgents,r_c)
% Hittar kollisioner genom hinder nu
% TODO hantera kollisioner på rätt sätt    
    
    col = zeros(3,1);
    for i = 1:numAgents-1
        for j = i+1:numAgents
            if norm(targets(i,:) - targets(j,:)) < r_c  

                line1 = [targets(i,:);pos(i,:)];
                line2 = [targets(j,:);pos(j,:)];
                
                collision_point = lineIntersection(line1,line2);

                col(1:2) = (targets(i,:)+targets(j,:))/2;
                col(3) = 1;
                rot_a(j) = 2*pi*rand;
                rot_a(i) = 2*pi*rand;
            end
        end
    end

end