
function [pos, rot_a, col] = moveAllAgents2(pos,targetPos, rot_a ,obstacle, L, threshold, r_c, mapSize)

    col = zeros(3,1);
    numAgents = size(pos,1);
    targets = zeros(numAgents,2);
    tempObstacle = obstacle;
    if ~isnan(mapSize)
        minPos = -floor(mapSize/2);
        cellIndex = floor((pos - minPos)/L) + 1;
    end
    
    for i = 1:numAgents
        %If in edge cell add walls to obstacle
       if ~isnan(mapSize)
           tempObstacle = obstacle;
           j = 1;
           if cellIndex(i,1) == 1 %v�nster
               tempObstacle(:,:,size(tempObstacle, 3) + j) = [0,0;0,1];
               j = j+1;
           end
           if cellIndex(i,1) == mapSize(1) %h�ger
               tempObstacle(:,:,size(tempObstacle, 3) + j) = [1,0;1,1];
               j = j+1;
           end
           if cellIndex(i,2) == 1 %ner
               tempObstacle(:,:,size(tempObstacle, 3) + j) = [0,0;1,0];
               j = j+1;
           end
           if cellIndex(i,2) == mapSize(2) %upp
               tempObstacle(:,:,size(tempObstacle, 3) + j) = [0,1;1,1];
           end
       end
       %Move agent
       [targets(i,:), rot_a(i)] = moveAgent(pos(i,:),targetPos(i,:), rot_a(i), tempObstacle, L, threshold);
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