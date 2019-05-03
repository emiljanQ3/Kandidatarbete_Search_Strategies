%%Iteration 2 of simulation
%CONFIG-------------------------------------------------------------------------------------------------------

obstacleType            = "c";
L                       = 30;      %Cell side length
mapSize                 = [3,3];  %Number of cells before a boundry is reached
edge                    = 1;  % true med kant fals utan
R                       = 0.6;
r                       = 0;
obstacle                = generateObstacle(obstacleType, R,r);   %Periodic obstacle contained in one cell
numAgents               = 1;
dT                      = 0.04;   % Delta time in seconds
preTime                 = 20;     %Number of seconds simulation is run before measurement starts.
measurmentTime          = 30;
numTimeSteps            = floor(measurmentTime/dT);
numSimulations          = 1000;
w                       = 10.^(linspace(-1,1,100));  % angle speed in rad/s      Should be defined as vector when doing tests for sevareal kiralities.
v                       = 16;     % 
l                       = 3.6; % Side length of cells in grid used to determine covered area
D_r                     = 0.04; %Diffusion constant for rotation
D_p                     = 0; %Diffusion constant for position
r_c                     = l;

%Config variables that might be interesting to include in the future:
%Friction
%Agent shape

%SETUP-----------------------------------------------------------------------------------------------------------------------

pos_a = zeros(numAgents, 2, numTimeSteps);      %INITIALIZATION: Agent positions in each timestep
pos_pre = zeros(numAgents, 2, floor(preTime/dT));
area = zeros(numSimulations,1);           %INITIALIZATION: List of the amount of area elements found each simulation.
totalTime = zeros(numSimulations,1); 
meanAreaCovered = zeros(length(w),1);           %INITIALIZATION: List of mean area covered for each kirality.
areaPerTime = zeros(length(w),1);               %INITIALIZATION: List of mean area covered per total time for each kirality.
colision = zeros(3,numTimeSteps);
tic
%SIMULATION LOOP-------------------------------------------------------------------------------------------------------------
w_j = 1;
w_count = 0;
loop_cycles = length(w);
startTic = tic;
for w_i = w %Loop over different kiralities
      
    for N_i = 1:numSimulations %Loop over separate simulations
        rot_a = 2*pi*rand(numAgents,1);          %Starting rotations
        pos_pre(:,:,1) = zeros(numAgents,2);     %Starting positions
        
        % Do the simulation for pre time
        [pos_pre, rot_a, colision, ~] = simulate( preTime ,dT,D_r,D_p,v,w_i,numAgents,pos_pre, rot_a,colision, obstacle, L, r_c,mapSize,false(1,1));

        %Set position to cell [1,1]
        pos_a(:,:,1) = pos_pre(:,:,end) - floor(pos_pre(:,:,end)/L)*L;
        
        % Do simulation for measurmenttime after pretime is done
        [pos_a, rot_a, colision,totalTime(N_i)] = simulate( measurmentTime ,dT,D_r,D_p,v,w_i,numAgents,pos_a, rot_a,colision, obstacle, L, r_c,mapSize,edge);
          
        %Simulation is done. Time to calculate area discovered.
        [area(N_i),~] = calcArea(pos_a,v,dT,l,1);
        
  
    end
    
    %All N simulations have been compleated. The mean result is saved for this kirality.
    meanAreaCovered(w_j) = mean(area);
    varians(w_j)        = std(area);
    areaPerTime(w_j) = mean(area./totalTime);
    w_j = w_j + 1;
    
    
    w_count = w_count + 1
        
    fprintf("Time elapsed:          " + sec2hms(toc(startTic)) + "\n")
    timeLeft = toc(startTic) / w_count * (loop_cycles - w_count);
    fprintf("Estimated time left:   " + sec2hms(timeLeft) + "\n\n")
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


%%
if(edge) 
    figure(112)
    areaPerTime_max = numAgents*(4*v*measurmentTime*l/pi-l^2)/measurmentTime;
    semilogx(w,areaPerTime/(areaPerTime_max),'o')
    axis([0.01, 10, 0, 1])
    name=strcat('step', num2str(dT), '; ', 'time', num2str(numTimeSteps), '; ', 'simulations', num2str(numSimulations), '; ', obstacleType, '; ', 'R=', num2str(R), '; D_r=', num2str(D_r));
    title(name)    
else
    figure(112)
    meanArea_max = numAgents*(4*v*measurmentTime*l/pi-l^2);
    semilogx(w,meanAreaCovered/(meanArea_max),'o')
    axis([0.01, 10, 0, 1.2])
    name=strcat('step', num2str(dT), '; ', 'time', num2str(numTimeSteps), '; ', 'simulations', num2str(numSimulations), '; ', obstacleType, '; ', 'R=', num2str(R), '; D_r=', num2str(D_r));
    title(name)
end

%% K�r detta script f�r att spara ditt workspace
dateTime = clock;
R_s = num2str(R);
r_s = num2str(r);
filename = strcat( join(string(dateTime(1:3)),''), '-', join(string(dateTime(4:5)),''), '_', obstacleType, R_s([1,3:end]), r_s([1,3:end]), '_', num2str(numAgents))
path = strcat(pwd, '/results/', filename)
save(path)
%saveas(h,figname, 'fig')

%% Animation of the last done kirality

p = animation(pos_a,obstacle,dT,colision, numTimeSteps);

%%

M = agentTracking('OG.xml')

r = splitPositionData(M)

%%

[cir, v] = getCirality(M,1/25,1)

[cir2, v2] = getKomplexCirality(r,1/25)

[cir3, v3] = getComplexCirality(r,1/25,1)





