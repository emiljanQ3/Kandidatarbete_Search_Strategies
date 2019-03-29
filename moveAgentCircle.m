function [resultPos, resultRotation] = moveAgentCircle(pos,targetPos, rotation, R)
    
    resultRotation = rotation;
    
    if(norm(targetPos) > R) %If moving outside circle, place on boundary
        resultPos = targetPos*R/norm(targetPos);
    else
        resultPos = targetPos;
    end
end
