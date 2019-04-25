function [resultPos, resultRotation] = moveAgentSquare(pos,targetPos, rotation, R)
    
    resultRotation = rotation;
    resultPos = targetPos;
    
    if (targetPos(1) > R) %If moving outside square, place on boundary
        resultPos(1) = R;
    elseif (targetPos(1) < -R)
        resultPos(1) = -R;
    end
    
    if (targetPos(2) > R) %If moving outside square, place on boundary
        resultPos(2) = R;
    elseif (targetPos(2) < -R)
        resultPos(2) = -R;
    end
end
