%% Reads a file, returns the calculated chirality and normalized area.

file = 'Filtered videos\2_Tracks.xml'; %Name of file
agent = 1; % which agent/s we look at

dT = 1/25; % Time step
M = cut(file,agent); % Turns file into a position matrix, without NaN:s

[kir,v] = getCirality(M,dT,1);

l = 30; % Corresponds to ~5*v*dT for agents in experiments
[squares,normA] = calcArea(M,v,dT,l);

result = [kir, normA]

%% Plot M, movement of agent

x = M(1,1,:);
y = M(1,2,:);
x = x(:,:)';
y = y(:,:)';
plot(x,y)

%% Run animation

% Generate an obstacle from the lab
obstacle = [0 0; 0 0];

p = animateExperiment(M,obstacle,dT);

%% Run to save workspace
dateTime = clock;
R_s = num2str(R);
r_s = num2str(r);
filename = strcat( join(string(dateTime(1:3)),''), '-', join(string(dateTime(4:5)),''), '_', obstacleType, R_s([1,3:end]), r_s([1,3:end]), '_', num2str(numAgents))
path = strcat(pwd, '/results/', filename)
save(path)
%saveas(h,figname, 'fig')

%% Plot graph for all chiralities
figure(112)
%areaPerTime_max = numAgents*(4*v*totalTime*l/pi-l^2)/totalTime;
areaPerTime_max = normA; %
semilogx(w,areaPerTime/(areaPerTime_max),'o')
axis([0.01, 10, 0, 1.2])
%name=strcat('step', num2str(dT), '; ', 'time', num2str(numTimeSteps), '; ', 'simulations', num2str(numSimulations), '; ', obstacleType, '; ', 'R=', num2str(R), '; D_r=', num2str(D_r));
name = 'Experimental results';
title(name)