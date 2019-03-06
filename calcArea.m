function [numSquares,A] = calcArea(M,v,dT,l)
% Calculates the area covered by an agent, returns the number of squares
% the agent moves through (numSquares) and the normalized area (A)
% Also consideres squares that are missed out when moving diagonally.
% The input Matrix should not contain any NaN, should be complete
% trajectories.

% Line 1:4 refers to the lines of a square, where line 1 is the left
% vertical line and the others are numbered in clockwise direction.

maxPos = max(max(M,[],1),[],3);
minPos = min(min(M,[],1),[],3);

gridSize = ceil((maxPos - minPos) / l);
areaGrid = zeros(gridSize);

indexedPos_a = floor((M - minPos)/l) + 1;

extraX = zeros(size(indexedPos_a,3),1);
extraY = zeros(size(indexedPos_a,3),1);

% Look for missing squares due to diagonal movement
fill = 1;
for a = 1:size(indexedPos_a,1)
    for i = 1:size(indexedPos_a,3)-1
       
       xMove = indexedPos_a(a,1,i)-indexedPos_a(a,1,i+1); 
       yMove = indexedPos_a(a,2,i)-indexedPos_a(a,2,i+1);

       if (norm(xMove)>0 && norm(yMove)>0)
           
           lineMovement = [M(a,1,i) M(a,2,i); M(a,1,i+1) M(a,2,i+1)];
           ind = indexedPos_a(a,:,i);
           
           lines = zeros(2,2,2);
           if (xMove > 0)
               % Generate line1
               lines(:,:,1) = [(ind(1)-1)*l+minPos(1) (ind(2)-1)*l+minPos(2); (ind(1)-1)*l+minPos(1) (ind(2))*l+minPos(2)];;
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
                  end
                  extraX(fill) = xInd;
                  extraY(fill) = yInd;
                  fill = fill + 1;
                  break;
                  
              end
           end
       end
    end
end

xIndices = indexedPos_a(:,1,:);
yIndices = indexedPos_a(:,2,:);
xIndices = xIndices(:,:)';
yIndices = yIndices(:,:)';

for i = 1:size(xIndices,1)
    areaGrid(xIndices(i),yIndices(i)) = 1;
    if (extraX(i) ~= 0 || extraY(i) ~= 0)
        areaGrid(extraX(i),extraY(i)) = 1;
    end
end

numSquares = sum(sum(areaGrid));
T = size(M,3)*dT;
% Normalized area
A = numSquares/((4*v*T/(pi*l)-1));


end