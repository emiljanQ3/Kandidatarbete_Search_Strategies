function [resultPos, resultRotation] = moveAgent2(pos,targetPos, rotation, obstacle,L, threshold)

    %If movement is small enough, return current position.
    if(norm(pos-targetPos) < threshold)
        resultPos = pos;
        resultRotation = rotation;
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
        resultRotation = rotation;
        return 
    end

    tangent = wall(2,:) - wall(1,:);
    tangent = tangent/norm(tangent);
    normal = tangent*[0,-1;1,0];

    
    %Calculate new target pos based on random collision direction.
    
    dotProdukt = (targetPos(1)-impactPoint(1))*normal(1) + (targetPos(2)-impactPoint(2))*normal(2);
    
    %==== USE ONE OF THE WAYS TO CALCULATE randomDirection ===
    
    %randomDirection = 2*pi*rand;
    %---------------------------
    %randomDirection = rotation + randn * 0.2;
    %---------------------------
    projectedPoint = pos + dot((impactPoint-pos), normal)*normal;
    rotSign = -sign((projectedPoint(2) - pos(2))*(impactPoint(1) - projectedPoint(1)) - (impactPoint(2) - projectedPoint(2))*(projectedPoint(1) - pos(1)));
    randomDirection = rotation + abs(randn) * 0.2 * rotSign;
    %---------------------------
%     v = (impactPoint - pos);
%     u = v - dot(v, normal)*2*normal;
%     dir = (u/norm(u) + tangent/norm(tangent) * sign(dot(u,tangent)))/2;
%     newRot = atan(dir(2)/dir(1));
%     if(dir(1) < 0)
%         newRot = newRot + pi;
%     end
%     randomDirection = newRot + randn * 0.2;
    %==========================================================
    
    lengthLeft = norm(pos-targetPos) - norm(impactPoint-targetPos);
    
    newTarget = impactPoint + [cos(randomDirection), sin(randomDirection)] * lengthLeft;
    newStartPos = impactPoint - 0.001*normal*dotProdukt;
    
    [resultPos, resultRotation] = moveAgent(newStartPos, newTarget, randomDirection, obstacle, L, threshold);
end
