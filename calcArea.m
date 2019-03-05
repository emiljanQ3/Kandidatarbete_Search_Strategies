function [numSquares,A] = calcArea(M,v,dT,l)

% Calculates the area covered by an agent, returns the number of squares
% the agent moves through (numSquares) and the normalized area (A)
% Also consideres squares that are missed out when moving diagonally.
% The input Matrix should not contain any NaN, should be complete
% trajectories.

%l = 30; % corresponds to ish 5*dT*v
dT_times_v = v*dT;

maxPos = max(max(M,[],1),[],3);
minPos = min(min(M,[],1),[],3);

gridSize = ceil((maxPos - minPos) / l);
areaGrid = zeros(gridSize);

%lostSquares = 0;
indexedPos_a = floor((M - minPos)/l) + 1;

extraX = zeros(size(indexedPos_a,3),1);
extraY = zeros(size(indexedPos_a,3),1);
fill = 1;
for a = 1:size(indexedPos_a,1)
    for i = 1:size(indexedPos_a,3)-1
       if (norm(indexedPos_a(a,1,i)-indexedPos_a(a,1,i+1))>0 && norm(indexedPos_a(a,2,i)-indexedPos_a(a,2,i+1))>0)
           lineMovement = [M(a,1,i) M(a,2,i); M(a,1,i+1) M(a,2,i+1)];
           ind = indexedPos_a(a,:,i);
           lines = zeros(2,2,1);
           lines(:,:,1) = [(ind(1)-1)*l+minPos(1) (ind(2)-1)*l+minPos(2); (ind(1)-1)*l+minPos(1) (ind(2))*l+minPos(2)];
           lines(:,:,2) = [(ind(1)-1)*l+minPos(1) (ind(2))*l+minPos(2); ind(1)*l+minPos(1) ind(2)*l+minPos(2)];
           lines(:,:,3) = [ind(1)*l+minPos(1) ind(2)*l+minPos(2); ind(1)*l+minPos(1) (ind(2)-1)*l+minPos(2)];
           lines(:,:,4) = [ind(1)*l+minPos(1) (ind(2)-1)*l+minPos(2); (ind(1)-1)*l+minPos(1) (ind(2)-1)*l+minPos(2)];
           
           for j = 1:4
               line = lines(:,:,j);
               %Reorder dimentions
               line = line(:,:);
              if (~isnan(lineIntersection(lineMovement,line)))
                  if (j == 1) 
                      xInd = ind(1)-1;
                      yInd = ind(2);
                  elseif j == 2
                      xInd = ind(1);
                      yInd = ind(2)+1;   
                  elseif j == 3
                      xInd = ind(1)+1;
                      yInd = ind(2);
                  elseif j == 4
                      xInd = ind(1);
                      yInd = ind(2)-1;
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
%lines(:,:,1);
xIndices = indexedPos_a(:,1,:);
yIndices = indexedPos_a(:,2,:);
xIndices = xIndices(:,:)';
yIndices = yIndices(:,:)';

for i = 1:size(xIndices,1)
    xIndices(i);
    yIndices(i);
    areaGrid(xIndices(i),yIndices(i)) = 1;
    i;
    if (extraX(i) ~= 0 || extraY(i) ~= 0)
        
        extraX(i);
        extraY(i);
        indexedPos_a(1,:,i:i+1);
        areaGrid(extraX(i),extraY(i)) = 1;
    end
end

numSquares = sum(sum(areaGrid)); 
%A = l*l*squares;
numTimeSteps = size(M,3);
numOfAgents = size(M,1);
% v*dT = 6 pixlar
A = numSquares/((4*v*dT/(pi*l)-1)); 
end