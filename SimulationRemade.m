%Iteration 2 of simulation

%CONFIG-------------------------------------------------------------------------------------------------------
L                       = 1;      %Cell side length
obstacle                = generateObstacle("hej", R,r);   %Periodic obstacle contained in one cell
numAgents               = 1;
numTimeSteps            = 300;
numSimulations          = 1;
dT                      = 0.1;   % Delta time in seconds
w                       = 0*pi/10;  % angle speed in rad/s      Should be defined as vector when doing tests for sevareal kiralities.
v                       = 1;     % speed in m/s
l                       = 2 * Delta_time * speed; % Side length of cells in grid used to determine covered area
D_r                     = 0.01; %Diffusion constant for rotation

%Config variables that might be interesting to include in the future:
%Friction
%Agent shape
%D_t   %Diffusion constant for movement
%Boundry of environment

%SETUP-----------------------------------------------------------------------------------------------------------------------

pos_a = zeros(2,numAgents, numTimeSteps);   %INITIALIZATION: Agent positions in each timestep
areaCovered = zeros(numSimulations);        %INITIALIZATION: List of the amount of area elements found each simulation.
meanAreaCovered = zeros(length(w));         %INITIALIZATION: List of mean area covered for each kirality.

%SIMULATION LOOP-------------------------------------------------------------------------------------------------------------
for w_i = w %Loop over different kiralities
    
    for N_i = 1:numSimulations %Loop over separate simulations
        
        rot_a = 2*pi*rand(numAgents); %Starting rotations
        pos_a(:,:,1) = null;          %Starting positions
        
        for T_i = 2:numTimeSteps
            rot_a = mod(rot_a + dT * w_i + sqrt(2 * D_r * dT) * randn, 2  * pi); %Update agent rotation for all agents
            
            for agent = 1:numAgents
                pos_a(:,agent,T_i) = moveAgent(pos_a(:,agent,T_i-1), rot_a(agent), v*dT, obstacle, L);    %Move agent
            end
            
        end
       
        %Simulation is done. Time to calculate area discovered.
        max = max(max(pos_a, 2),3);
        min = min(min(pos_a, 2),3);
        size = (max - min) / l;
        areaGrid = zeros(size(1), size(2));     %Giovanni told us to use sparse grid but I'm not sure how it increases efficency.
        
        indexedPos_a = ceil((pos_a + min)/l);
        
        areaGrid(indexedPos(1,:,:), indexedPos(2,:,:)) = 1;     %Not sure about the dimensions here, will check back later.
        
        areaCovered(simulation) = sum(areaGrid, 'all');

    end
    
    %All N simulations have been compleated. The mean result is saved for this kirality.
    meanAreaCovered(w_i) = mean(areaCovered);
    
end

%Result is stored as data points, pairing each kirality with a meanAreaCovered value.
result(:,1) = w;
result(:,2) = meanAreaCovered;