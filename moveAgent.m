function resultPos = moveAgent(pos,targetPos,obstacle,L, threshold)

    %If movement is small enough, return current position.
    if(norm(pos-targetPos) < threshold)
        resultPos = pos;
        return
    end
    
    impactDist = inf;       %Distance to closest wall collision.
    impactPoint = NaN;      %Point of collision
    wall = zeros(2,2);      %Closest wall that causes a collision.

    for i = 1:size(obstacle,3)                              %TODO fix collision when entering new cell
        intersectPoint = lineIntersection([pos; targetPos] - floor(pos/L)*L, obstacle(:,:,i)*L); %Get intersection between wall segment and movement, in cell 1.
        if ~isnan(intersectPoint) %If intersection found
            intersectPoint = floor(pos/L)*L + intersectPoint; %Move impact point to correct cell.
            intersectDist = norm(pos-intersectPoint); 
            if intersectDist < impactDist %If intersection is closer than closest discovered impact
                %Set intersection as new impact point
                impactDist = intersectDist; 
                impactPoint = intersectPoint;
                wall = obstacle(:,:,i);
            end
        end
    end

    if(isnan(impactPoint)) %If no wall impact, move unhindered.
        resultPos = targetPos; 
        return 
    end

    tangent = wall(2,:) - wall(1,:);
    tangent = tangent/norm(tangent);
    normal = tangent*[0,-1;1,0];

    
    %Calculate new target pos based on tangental movement along the wall.
    
    dotProdukt = (targetPos(1)-impactPoint(1))*normal(1) + (targetPos(2)-impactPoint(2))*normal(2);
    newTarget = targetPos - 1.001*normal*dotProdukt;
    newStartPos = impactPoint - 0.001*normal*dotProdukt;
    
    % 
    resultPos = moveAgent(newStartPos, newTarget, obstacle, L, threshold);
end
