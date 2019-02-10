function resultPos = moveAgent(pos,targetPos,obstacle,L, threshold)
    
    impactDist = inf;       %Distance to closest wall collision.
    impactPoint = NaN;      %Point of collision
    wall = zeros(2,2);      %Closest wall that causes a collision.

    for i = 1:size(obstacle,3)                              %TODO fix collision when entering new cell
        intersectPoint = lineIntersection([pos; targetPos] - floor(pos/L)*L, obstacle(:,:,i)*L); %Get intersection between wall segment and movement, in cell 1.
        if ~isnan(intersectPoint) %If intersection found
            intersectPoint = floor(pos/L)*L + intersectPoint; %Move impact point to correct cell.
            %--DEBUG----------
            scatter(intersectPoint(1), intersectPoint(2), 'r')
            scatter(pos(1),pos(2), 'y')
            resultPos = pos;
            return
            %--END DEBUG-----------
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
    targetPos = impactPoint + tangent*dot((pos-impactPoint), tangent) - 0.00001*normal*dot((pos-impactPoint), normal);

    %If movement is small enough, return current position.
    if(norm(impactPoint-targetPos) < threshold)
        resultPos = pos;
        return
    end

    resultPos = moveAgent(impactPoint, targetPos, obstacle, L, threshold);
end
