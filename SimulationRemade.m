%%Iteration 2 of simulation
clf
close all
hold on
warning('off','all')
%CONFIG-------------------------------------------------------------------------------------------------------
obstacleType            = "c";
L                       = 1;      %Cell side length
R                       = 0.5;
r                       = 0.05;
obstacle                = generateObstacle(obstacleType, R,r);   %Periodic obstacle contained in one cell
numAgents               = 1;
numTimeSteps            = 300;
numSimulations          = 10;

dT                      = 0.1;   % Delta time in seconds
w                       = 10.^linspace(-2,1,100);  % angle speed in rad/s      Should be defined as vector when doing tests for sevareal kiralities.
v                       = 1;     % speed in m/s
l                       = 1.5 * dT * v; % Side length of cells in grid used to determine covered area
D_r                     = 0.01; %Diffusion constant for rotation

%Config variables that might be interesting to include in the future:
%Friction
%Agent shape
%D_t   %Diffusion constant for movement
%Boundry of environment

%SETUP-----------------------------------------------------------------------------------------------------------------------

pos_a = zeros(numAgents, 2, numTimeSteps);  %INITIALIZATION: Agent positions in each timestep
areaCovered = zeros(numSimulations,1);        %INITIALIZATION: List of the amount of area elements found each simulation.
<<<<<<< HEAD
meanAreaCovered = zeros(length(w));         %INITIALIZATION: List of mean area covered for each kirality.
=======
meanAreaCovered = zeros(length(w),1);         %INITIALIZATION: List of mean area covered for each kirality.
tic
%SIMULATION LOOP-------------------------------------------------------------------------------------------------------------
w_j = 1;
for w_i = w %Loop over different kiralities
    w_j
    for N_i = 1:numSimulations %Loop over separate simulations
        
        rot_a = 2*pi*rand(numAgents); %Starting rotations
        pos_a(:,:,1) = zeros;          %Starting positions
        
        for T_i = 2:numTimeSteps
            rot_a = mod(rot_a + dT * w_i + sqrt(2 * D_r * dT) * randn, 2  * pi); %Update agent rotation for all agents
            
            for agent = 1:numAgents
                targetPos = pos_a(agent, :, T_i-1) + [cos(rot_a(agent)), sin(rot_a(agent))] * dT * v; %Calculate where a unhindered move would go.
                pos_a(agent, :, T_i) = moveAgent(pos_a(agent, :, T_i-1), targetPos, obstacle, L, v*dT/10);    %Move agent and take obstacles into consideration.
            end
            
        end 
        
        %Simulation is done. Time to calculate area discovered.
        maxPos = max(max(pos_a,[],1),[],3);
        minPos = min(min(pos_a,[],1),[],3);
        gridSize = ceil((maxPos - minPos) / l);
        areaGrid = zeros(gridSize);
         

         indexedPos_a = ceil((pos_a - minPos)/l)+1; %so as to start at index 1
        %Reorder dimentions
         for i = 1:size(indexedPos_a,3)                %HELP lyckas inte med vektorer som index 
             areaGrid(indexedPos_a(agent,1,i),indexedPos_a(agent,2,i)) = 1;
        end   
 %        areaGrid
       % areaGrid(indexedPos_a(agent, 1, :),indexedPos_a(agent, 2, :)) = 1;     %Not sure about the dimensions here, will check back later.         
       areaCovered(N_i) = sum(sum(areaGrid));

        indexedPos_a = floor((pos_a - minPos)/l) + 1;
        
        xIndices = indexedPos_a(:,1,:);
        yIndices = indexedPos_a(:,2,:);
        for i = 1:numel(xIndices)             
            areaGrid(xIndices(i),yIndices(i)) = 1;
        end   
              
        areaCovered(N_i) = sum(areaGrid,'all');
       

    end
    
    %All N simulations have been compleated. The mean result is saved for this kirality.
<<<<<<< HEAD
    meanAreaCovered(1) = mean(areaCovered);
=======
    meanAreaCovered(w_j) = mean(areaCovered);
    w_j = w_j + 1;
>>>>>>> 20b1ea59f6b689888c2b295b079005369277280c
    
end

%Result is stored as data points, pairing each kirality with a meanAreaCovered value.
semilogx(w,meanAreaCovered)
%Plot-----------------------------------------------------------------------------------------------------------------
figure(101)
hold on
plotSize = 3;
for i = -plotSize:plotSize
    for j = -plotSize:plotSize
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
toc

%% KÖr detta script för att spara ditt workspace
dateTime = clock;
R_s = num2str(R);
r_s = num2str(r);
filename = strcat( join(string(dateTime(1:3)),''), '-', join(string(dateTime(4:5)),''), '_', obstacleType, R_s([1,3:end]), r_s([1,3:end]), '_', num2str(numAgents))
path = strcat(pwd, '\results\', filename)
save(path)
