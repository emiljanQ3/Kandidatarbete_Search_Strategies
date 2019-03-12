%% Reads a file, returns the calculated chirality and normalized area.

%file = 'XMLfiles/Kors_1agent/011_Tracks.xml'; %Name of file
agent = 1; % which agent/s we look at
l = 15; % Corresponds to ~5*v*dT for agents in experiments

dT = 1/25; % Time step

file = 'XMLfiles/Kors_1agent/013_Tracks.xml';

[pos_a,length,times] = cut(file,agent); % Turns file into a position matrix, without NaN:s

r = splitPositionData(pos_a)

[kir,v] = getComplexCirality(r,dT,1);

[squares,normA] = calcArea(pos_a,v,dT,l);

result = [kir, normA]

v

%% Run to save results to file
clc;
% Sparar:
    % Filnamn/path (i en egen fil)
    % Kiralitet
    % NormArea/tid
    % Total tid
    % hastighet
    % l (size of area elements)

expName = 'c1agent'; %Change name for each new set of data

file1 = ['results/Lab/' expName '.txt']; % Name of dataFile
file2 = ['results/Lab/' expName 'SourceFiles.txt']; % Name of file containing names of XML files
    
data = [result size(pos_a,3)*dT v l];
dlmwrite(file1,data,'-append');

fileID = fopen(file2,'a');
fprintf(fileID,'%-40s\n',file);
fclose(fileID);

%% Plot M, trajectory of agent

x = pos_a(1,1,:);
y = pos_a(1,2,:);
x = x(:,:)';
y = y(:,:)';
plot(x,y)
axis equal
%% Run animation

% Generate an obstacle from the lab
obstacle = [0 0; 0 0];

p = animateExperiment(pos_a,obstacle,dT);


%% Plot graph for all chiralities

experiment = 'hm1agent';
dataFile = ['results/Lab/' experiment '.txt'];
sourceFile = ['results/Lab/' experiment 'SourceFiles.txt'];
allData = dlmread(dataFile);

w = abs(allData(:,1));
normArea = allData(:,2); 
v = allData(:,4);



figure(113)
semilogx(w,normArea,'o')
figure(114)
semilogx(w./v,normArea,'o')

%axis([0.01, 10, 0, 1.2])
%name=strcat('step', num2str(dT), '; ', 'time', num2str(numTimeSteps), '; ', 'simulations', num2str(numSimulations), '; ', obstacleType, '; ', 'R=', num2str(R), '; D_r=', num2str(D_r));
name = 'Experimental results';
title(name)






