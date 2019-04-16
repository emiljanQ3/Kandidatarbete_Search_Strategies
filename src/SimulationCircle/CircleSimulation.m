%%Iteration 2 of simulation
%CONFIG-------------------------------------------------------------------------------------------------------

R                       = 1;      %Circular areana radius
numAgents               = 1;
dT                      = 0.04;   % Delta time in seconds
preTime                 = 5;     % Number of seconds simulation is run before measurement starts.
measurmentTime          = 50;
numTimeSteps            = floor(measurmentTime/dT);
numSimulations          = 1000;
w                       = 10.^(linspace(-2,1,150));  % angle speed in rad/s      Should be defined as vector when doing tests for sevareal kiralities.
v                       = 0.652;     % speed in R/s
l                       = 0.156; % Side length of cells in grid used to determine covered area
D_r                     = 0.02; %Diffusion constant for rotation
D_p                     = 0; %Diffusion constant for position
r_c                     = l/2;
numAreaDP               = 100;

%Config variables that might be interesting to include in the future:
%Friction
%Agent shape
%SETUP-----------------------------------------------------------------------------------------------------------------------

pos_a = zeros(numAgents, 2, numTimeSteps);              %INITIALIZATION: Agent positions in each timestep
pos_pre = zeros(numAgents, 2, floor(preTime/dT));
area = zeros(numSimulations,numAreaDP);           %INITIALIZATION: List of the amount of area elements found each simulation.
normA      = zeros(numSimulations,numAreaDP);
totalTime = zeros(numSimulations,1); 
meanAreaCovered = zeros(length(w),numAreaDP);           %INITIALIZATION: List of mean area covered for each kirality.
areaPerTime = zeros(length(w),numAreaDP);               %INITIALIZATION: List of mean area covered per total time for each kirality.
colision = zeros(3,numTimeSteps);
tic
%SIMULATION LOOP-------------------------------------------------------------------------------------------------------------
w_j = 1;
for w_i = w %Loop over different kiralities
    w_j
      
    for N_i = 1:numSimulations %Loop over separate simulations
        rot_a = 2*pi*rand(numAgents,1);          %Starting rotations
        pos_pre(:,:,1) = zeros(numAgents,2);     %Starting positions
        
        % Do the simulation for pre time
        [pos_pre, rot_a, colision, ~] = simulateCircle( preTime ,dT,D_r,D_p,v,w_i,numAgents,pos_pre, rot_a,colision, R, r_c);
        
        %Set starting pos for measurement
        pos_a(:,:,1) = pos_pre(:,:,end);
        
        % Do simulation for measurmenttime after pretime is done
        [pos_a, rot_a, colision,totalTime(N_i)] = simulateCircle( measurmentTime ,dT,D_r,D_p,v,w_i,numAgents,pos_a, rot_a,colision, R, r_c);
        
        %Let's also calculate area discovered at certain times

        [area(N_i, :), normA(N_i, :)] = calcArea(pos_a, v, dT, l, numAreaDP);

    end
    
    %All N simulations have been compleated. The mean result is saved for this kirality.
    meanNormA(w_j,:)       = mean(normA);
    meanAreaCovered(w_j,:) = mean(area);
    varians(w_j)        = std(area(:, numAreaDP));
    %areaPerTime(w_j,:) = meanAreaCovered(w_j,:)./totalTime;
    w_j = w_j + 1;
    toc
end
toc

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
l_s = num2str(l);
time_s = num2str(measurmentTime);

filename = strcat( join(string(dateTime(1:3)),''), '-', join(string(dateTime(4:5)),''), '_', 'circle', '_R', R_s([1,3:end]), '_t', time_s, '_l', l_s([1,3:end]));
path = strcat(pwd, '/results/Final_simulation/', filename)
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





