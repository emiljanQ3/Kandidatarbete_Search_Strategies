function [numSquares,A] = calcArea(M,vdT)

% Calculates the area covered by an agent, returns the number of squares
% the agent moves through (numSquares) and the normalized area (A) 
% The input Matrix should not contain any NaN, should be complete
% trajectories.

% dT = 0.04; % 1/25
% v = 0;
% for i = 1:size(M,1)
%     [k, v_ny] = getCirality(M(i,:,2:11),t_step);
%     if v_ny > v
%         v = v_ny
%     end
% end
%l = 1.5 * dT * v;

l = 30; % corresponds to 5*dT*v
dT_times_v = vdT;

maxPos = max(max(M,[],1),[],3);
minPos = min(min(M,[],1),[],3);

gridSize = ceil((maxPos - minPos) / l);
areaGrid = zeros(gridSize);

indexedPos_a = floor((M - minPos)/l) + 1;
indexedNoFloor = (M - minPos)/l;

xIndices = indexedPos_a(:,1,:);
yIndices = indexedPos_a(:,2,:);
xIndices = xIndices(:,:)';
yIndices = yIndices(:,:)';

for i = 1:size(xIndices)
    xIndices(i);
    yIndices(i);
    areaGrid(xIndices(i),yIndices(i)) = 1;
end

numSquares = sum(sum(areaGrid)); 
%A = l*l*squares;
numTimeSteps = size(M,3)
numOfAgents = size(M,1);
% v*dT = 6 pixlar
A = numSquares/(numOfAgents*dT_times_v*numTimeSteps/l)
end