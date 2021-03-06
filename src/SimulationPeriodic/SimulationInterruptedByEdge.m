%%Iteration 2 of simulation
%CONFIG-------------------------------------------------------------------------------------------------------

obstacleType            = "c";
L                       = 1;      %Cell side length
mapSize                 = [3,3]; %Number of cells before a boundry is reached
R                       = 0.8;
r                       = 0.5;
obstacle                = generateObstacle(obstacleType, R,r);   %Periodic obstacle contained in one cell
numAgents               = 1;
dT                      = 0.1;   % Delta time in seconds
preTime                 = 10;     %Number of seconds simulation is run before measurement starts.
simulationTime          = 15;
numTimeSteps            = floor(simulationTime/dT);
numSimulations          = 1000;
w                       = 10.^(linspace(-1,1,100));  % angle speed in rad/s      Should be defined as vector when doing tests for sevareal kiralities.
v                       = 1;     % speed in m/s
l                       = 0.156; % Side length of cells in grid used to determine covered area
D_r                     = 0.01; %Diffusion constant for rotation
D_p                     = 0; %Diffusion constant for position
r_c                     = l/2;

%Config variables that might be interesting to include in the future:
%Friction
%Agent shape

%SETUP-----------------------------------------------------------------------------------------------------------------------

pos_a = zeros(numAgents, 2, numTimeSteps);    %INITIALIZATION: Agent positions in each timestep
pos_pre = zeros(numAgents, 2, floor(preTime/dT));
areaCovered = zeros(numSimulations,1);        %INITIALIZATION: List of the amount of area elements found each simulation.
meanAreaCovered = zeros(length(w),1);         %INITIALIZATION: List of mean area covered for each kirality.
colision = zeros(3,numTimeSteps);
tic
%SIMULATION LOOP-------------------------------------------------------------------------------------------------------------
w_j = 1;
for w_i = w %Loop over different kiralities
      w_j
    for N_i = 1:numSimulations %Loop over separate simulations
        
        rot_a = 2*pi*rand(numAgents,1); %Starting rotations
        pos_pre(:,:,1) = zeros(numAgents,2); %randn(numAgents,2);          %Starting positions
        
        %Pre measurement start
        for t_i = 2:floor(preTime/dT)
            rot_a = mod(rot_a + dT * w_i + sqrt(2 * D_r * dT) * randn(size(rot_a)), 2  * pi); %Update agent rotation for all agents
            targetPos = pos_pre(:, :, t_i-1) + [cos(rot_a), sin(rot_a)] * dT * v + randn(numAgents, 2) * sqrt(2 * D_p * dT); %Calculate where a unhindered move would go.
            [pos_pre(:, :, t_i), rot_a, col]= moveAllAgents(pos_pre(:, :, t_i-1), targetPos,rot_a, obstacle, L, v*dT/10, r_c, mapSize);    %Move agent and take obstacles into consideration.        
        end
        
        %Set position to cell [1,1]
        pos_a(:,:,1) = pos_pre(:,:,end) - floor(pos_pre(:,:,end)/L)*L;
        
        %Simulations for measured values
        for T_i = 2:numTimeSteps
            rot_a = mod(rot_a + dT * w_i + sqrt(2 * D_r * dT) * randn(size(rot_a)), 2  * pi); %Update agent rotation for all agents
            targetPos = pos_a(:, :, T_i-1) + [cos(rot_a), sin(rot_a)] * dT * v + randn(numAgents, 2) * sqrt(2 * D_p * dT); %Calculate where a unhindered move would go.
            [pos_a(:, :, T_i), rot_a, col]= moveAllAgents(pos_a(:, :, T_i-1), targetPos,rot_a, obstacle, L, v*dT/10, r_c, mapSize);    %Move agent and take obstacles into consideration.        
            colision(:,T_i) = col;
            
            %Check if outside
            minPos = -floor(mapSize/2);
            cellIndex = floor((pos_a(:,:,T_i) - minPos)/L) + 1;
            if any(cellIndex(1) > mapSize(1) || cellIndex(2) > mapSize(2) || cellIndex(1) < 1 || cellIndex(2) < 1)
               break
            end
        end 
        
        %Simulation is done. Time to calculate area discovered.
        maxPos = max(max(pos_a,[],1),[],3);
        minPos = min(min(pos_a,[],1),[],3);
        gridSize = ceil((maxPos - minPos) / l);
        areaGrid = zeros(gridSize);
         
        indexedPos_a = floor((pos_a - minPos)/l) + 1;
        
        xIndices = indexedPos_a(:,1,:);
        yIndices = indexedPos_a(:,2,:);
        for i = 1:numel(xIndices)             
            areaGrid(xIndices(i),yIndices(i)) = 1;
        end   
              
        areaCovered(N_i) = sum(sum(areaGrid));
  
    end
    
    %All N simulations have been compleated. The mean result is saved for this kirality.
    meanAreaCovered(w_j) = mean(areaCovered);
    w_j = w_j + 1;
    
end
toc
%Result is stored as data points, pairing each kirality with a meanAreaCovered value.
%semilogx(w,meanAreaCovered)

%%
%Plot-----------------------------------------------------------------------------------------------------------------
hold on

maxPos = max(max(pos_a,[],1),[],3);
minPos = min(min(pos_a,[],1),[],3);
maxmax = max(abs([maxPos minPos]));
plotSize = ceil((maxmax) / L);
for i = -plotSize():plotSize-1
    for j = -plotSize:plotSize-1
        for k = 1:size(obstacle, 3)
            plot(obstacle(:,1,k)+j*L,obstacle(:,2,k)+i*L, 'k', 'LineWidth', 1)
            
        end
    end
end

for agent = 1:numAgents
    X = pos_a(agent, 1, :);
    Y = pos_a(agent, 2, :);
    %Reorder dimentions
    X = X(:,:)';
    Y = Y(:,:)';
    plot(X,Y);
end
%scatter(pos_a(agent, 1, :),pos_a(agent, 2, :), 'b.')

%%
figure(112)
semilogx(w,meanAreaCovered/(numAgents*v*dT*numTimeSteps/l),'o')
axis([0.01, 10, 0, 1.2])
name=strcat('step', num2str(dT), '; ', 'time', num2str(numTimeSteps), '; ', 'simulations', num2str(numSimulations), '; ', obstacleType, '; ', 'R=', num2str(R), '; D_r=', num2str(D_r));
title(name)
%% K�r detta script f�r att spara ditt workspace
dateTime = clock;
R_s = num2str(R);
r_s = num2str(r);
filename = strcat( join(string(dateTime(1:3)),''), '-', join(string(dateTime(4:5)),''), '_', obstacleType, R_s([1,3:end]), r_s([1,3:end]), '_', num2str(numAgents))
path = strcat(pwd, '/results/', filename)
save(path)
%saveas(h,figname, 'fig')

%% Animation of the last done kirality

p = animation(pos_a,obstacle,dT,colision);


