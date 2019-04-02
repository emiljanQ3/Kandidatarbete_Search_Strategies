function [numSquares,A] = calcArea(M,v,dT,l,nReturn)
% Calculates the area covered by an agent, returns the number of squares
% the agent moves through (numSquares) and the normalized area (A)
% Also consideres squares that are missed out when moving diagonally.
% The input Matrix should not contain any NaN, should be complete
% trajectories.

% Line 1:4 refers to the lines of a square, where line 1 is the left
% vertical line and the others are numbered in clockwise direction.

numSquares = zeros(1,nReturn);
A = zeros(1,nReturn);

maxPos = max(max(M,[],1),[],3);
minPos = min(min(M,[],1),[],3);

gridSize = ceil((maxPos - minPos) / l);
areaGrid = zeros(gridSize);

indexedPos_a = floor((M - minPos)/l) + 1;

extraX = zeros(size(indexedPos_a,3),1);
extraY = zeros(size(indexedPos_a,3),1);

% Look for missing squares due to diagonal movement
for a = 1:size(indexedPos_a,1)
    for i = 1:size(indexedPos_a,3)-1
        
        xMove = indexedPos_a(a,1,i)-indexedPos_a(a,1,i+1);
        yMove = indexedPos_a(a,2,i)-indexedPos_a(a,2,i+1);
        
        if ~(norm(xMove)>0 && norm(yMove)>0)
            continue;
        end
            
        lineMovement = [M(a,1,i) M(a,2,i); M(a,1,i+1) M(a,2,i+1)];
        ind = indexedPos_a(a,:,i);

        lines = zeros(2,2,2);
        if (xMove > 0)
            % Generate line1
            lines(:,:,1) = [(ind(1)-1)*l+minPos(1) (ind(2)-1)*l+minPos(2); (ind(1)-1)*l+minPos(1) (ind(2))*l+minPos(2)];
        else
            % Generate line3
            lines(:,:,1) = [ind(1)*l+minPos(1) ind(2)*l+minPos(2); ind(1)*l+minPos(1) (ind(2)-1)*l+minPos(2)];
        end

        if (yMove > 0)
            % Generate line4
            lines(:,:,2) = [ind(1)*l+minPos(1) (ind(2)-1)*l+minPos(2); (ind(1)-1)*l+minPos(1) (ind(2)-1)*l+minPos(2)];
        else
            % Generate line2
            lines(:,:,2) = [(ind(1)-1)*l+minPos(1) (ind(2))*l+minPos(2); ind(1)*l+minPos(1) ind(2)*l+minPos(2)];
        end

        % Check if agent crosses horizontal or vertical line
        for j = 1:2
            line = lines(:,:,j);
            %Reorder dimentions
            %line = line(:,:);
            if (~isnan(lineIntersection(lineMovement,line)))
                if (j == 1) % Crosses vertical lines
                    xInd = ind(1)-xMove; % xMove negative if positive movement
                    yInd = ind(2);
                elseif j == 2 % Crosses horizontal lines
                    xInd = ind(1);
                    yInd = ind(2)-yMove;
                elseif j == 3
                    xInd = ind(1)+1;
                    yInd = ind(2);
                elseif j == 4
                    xInd = ind(1);
                    yInd = ind(2)-1;
                end

                extraX(i) = xInd;
                extraY(i) = yInd;
                break;

            end
        end
    end
end

xIndices = indexedPos_a(:,1,:);
yIndices = indexedPos_a(:,2,:);
xIndices = xIndices(:,:)';
yIndices = yIndices(:,:)';

returnIndices = floor(linspace(1,size(xIndices,1),nReturn));
k = 1;
for j = 1:nReturn
    
    for i = k:returnIndices(j)
        
        areaGrid(xIndices(i),yIndices(i)) = 1;
        
        if (extraX(i) ~= 0 || extraY(i) ~= 0)
            areaGrid(extraX(i),extraY(i)) = 1;
        end
        
    end
    
    numSquares(j) = sum(sum(areaGrid));
    
    T = returnIndices(j)*dT;
    % Normalized area
    A(j) = numSquares(j)/((4*v*T/(pi*l)-1));
    
    k = returnIndices(j);
end



end