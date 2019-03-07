%% Reads a file, returns the calculated chirality and normalized area.

file = 'Filtered videos\2_Tracks.xml'; %Name of file
agent = 1; % which agent/s we look at

dT = 1/25; % Time step
[M,length,times] = cut(file,agent); % Turns file into a position matrix, without NaN:s

[kir,v] = getCirality(M,dT,1);

l = 30; % Corresponds to ~5*v*dT for agents in experiments
[squares,normA] = calcArea(M,v,dT,l);

result = [kir, normA]

%% Run to save results to file
clc;
% Sparar:
    % Filnamn/path (i en egen fil)
    % Kiralitet
    % NormArea
    % Total tid
    % hastighet
    % l (size of area elements)

expName = 'hm1agent'; %Change name for each new set of data

file1 = ['results\Lab\' expName '.txt']; % Name of dataFile
file2 = ['results\Lab\' expName 'SourceFiles.txt']; % Name of file containing names of XML files
    
data = [result size(M,3)*dT v l];
dlmwrite(file1,data,'-append');

fileID = fopen(file2,'a');
fprintf(fileID,'%-40s\n',file);
fclose(fileID);

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


%% Plot graph for all chiralities

experiment = 'hm1agent';
dataFile = ['results\Lab\' experiment '.txt'];
sourceFile = ['results\Lab\' experiment 'SourceFiles.txt'];
allData = dlmread(dataFile);

w = norm(allData(:,1));
normArea = allData(:,2); 

figure(112)
semilogx(w,normArea,'o')
%axis([0.01, 10, 0, 1.2])
%name=strcat('step', num2str(dT), '; ', 'time', num2str(numTimeSteps), '; ', 'simulations', num2str(numSimulations), '; ', obstacleType, '; ', 'R=', num2str(R), '; D_r=', num2str(D_r));
name = 'Experimental results';
title(name)