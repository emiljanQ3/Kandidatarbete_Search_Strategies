%%Iteration 2 of simulation
%CONFIG-------------------------------------------------------------------------------------------------------

R                       = 1;      %Circular areana radius
r                       = 0.167;
numAgents               = 2;
dT                      = 0.04;   % Delta time in seconds
preTime                 = 10;     %Number of seconds simulation is run before measurement starts.
maxMeasurmentTime       = 300;
numTimeSteps            = floor(maxMeasurmentTime/dT);
numSimulations          = 100;
w                       = linspace(-4,4,50); %[0.3366, 0.7897, 1.1479, 1.7525, 3.8640]; %10.^(linspace(-1,1,100));  % angle speed in rad/s      Should be defined as vector when doing tests for sevareal kiralities.
v                       = 0.5;     % speed in R/s
L                       = 1/7.5*R; % Side length of cells in grid used to determine covered area
D_r                     = 0.03; %Diffusion constant for rotation
D_p                     = 0; %Diffusion constant for position
r_c                     = R/5;
numAreaDP               = 100;

%SETUP-----------------------------------------------------------------------------------------------------------------------

pos_a = zeros(numAgents, 2, numTimeSteps);              %INITIALIZATION: Agent positions in each timestep
pos_pre = zeros(numAgents, 2, floor(preTime/dT));
numSquares = zeros(numSimulations,numAreaDP);           %INITIALIZATION: List of the amount of area elements found each simulation.
normA      = zeros(numSimulations,numAreaDP);
totalTimeSteps = zeros(numSimulations,1); 
meanAreaCovered = zeros(length(w),numAreaDP);           %INITIALIZATION: List of mean area covered for each kirality.
areaPerTime = zeros(length(w),numAreaDP);               %INITIALIZATION: List of mean area covered per total time for each kirality.
colision = zeros(3,numTimeSteps);


%SIMULATION LOOP-------------------------------------------------------------------------------------------------------------
w_count = 1;
loop_cycles = length(w)*(length(w)-1)/2;

startTic = tic;

for i = 1:length(w)
    for j = i:1:length(w)
        loopTic = tic;

        W = [w(i),w(j)];
        
        for N_i = 1:numSimulations %Loop over separate simulations

            rot_a = [pi; 0];
            pos_a(:,:,1) =  [-1/3, 0;
                              1/3, 0];

            % Do simulation for measurmenttime after pretime is done
            [pos_a, rot_a, colision, totalTimeSteps(N_i)] = simulateMultiAgentCircle( maxMeasurmentTime ,dT,D_r,D_p,v,W,numAgents,pos_a, rot_a,colision, R, r_c,1);

        end

        meanTotalTime(i,j) = mean(totalTimeSteps)*dT;


        status = string(w_count) + "/" + string(loop_cycles) + "   Chirality: [" + string(w(i)) + ", " + string(w(j)) + "]   Mean time: " + string(meanTotalTime(i,j))

        w_count = w_count + 1;
        toc(loopTic)
        %animation(pos_a, generateObstacle('hm'),0,colision, totalTime(N_i))
    end
end
toc(startTic)
%% 3D plot
X = w;
Y = w;
Z = meanTotalTime + meanTotalTime' - diag(diag(meanTotalTime));
Z_1 = Z;
Z_2 = 1./Z;

figure
surf(X,Y,Z_1)
title("Time")

figure
surf(X,Y,Z_2)
title("Efficiency")
%%
figure
hold on
c = jet(length(w));
for i = 1:size(meanAreaCovered,1)
    plot(1:numAreaDP, meanAreaCovered(i,:)./pi, 'color', c(i,:))
end

figure
hold on
for i = 1:size(meanNormA,1)
    plot(1:numAreaDP, meanNormA(i,:)./pi)
end

%%
%Plot-----------------------------------------------------------------------------------------------------------------
hold on
obstacle = generateObstacle('hm');

maxPos = max(max(pos_a,[],1),[],3);
minPos = min(min(pos_a,[],1),[],3);
maxmax = max(abs([maxPos minPos]));
plotSize = ceil((maxmax) / R);
for i = -plotSize():plotSize-1
    for j = -plotSize:plotSize-1
        for k = 1:size(obstacle, 3)
            plot(obstacle(:,1,k)+j*R,obstacle(:,2,k)+i*R, 'k', 'LineWidth', 1)
            
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



%% K�r detta script f�r att spara ditt workspace
dateTime = clock;
R_s = num2str(R);
r_s = num2str(r);
filename = strcat( join(string(dateTime(1:3)),''), '-', join(string(dateTime(4:5)),''), '_', 'circle', R_s([1,3:end]), r_s([1,3:end]), '_', num2str(numAgents))
path = strcat(pwd, '/results/', filename)
save(path)
%saveas(h,figname, 'fig')

%% Animation of the last done kirality

p = animation(pos_a, obstacle,dT,colision);

%%

M = agentTracking('OG.xml')

r = splitPositionData(M)

%%

[cir, v] = getCirality(M,1/25,1)

[cir2, v2] = getKomplexCirality(r,1/25)

[cir3, v3] = getComplexCirality(r,1/25,1)





