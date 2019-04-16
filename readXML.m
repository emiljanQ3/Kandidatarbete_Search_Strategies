%% Reads a file, returns the calculated chirality and normalized area.

file = 'XMLfiles/Circle_large_1agent/Stor (54)_Tracks.xml'; %Name of file

agent = 1; % which agent/s we look at

l = 10; % Corresponds to ~5*v*dT for agents in experiments

dT = 1/25; % Time step

[pos_a,length,times] = cut(file,agent); % Turns file into a position matrix, without NaN:s

[r, indice] = splitPositionDataPartitioned(pos_a,700, myCircle);

[kir,D_r,v] = getKompSpiral(r,dT,1,6,60);
indice
%[squares,normA] = calcArea(pos_a,v,dT,l,size(pos_a,3),1);

result = [kir D_r v size(pos_a,3)*dT]%, normA]

%% Run to get data from circular path (run once)
file = 'XMLfiles/Circle_large_1agent/Stor (54)_Tracks.xml'; %Name of file

agent = 1; % which agent/s we look at

l = 10; % Corresponds to ~5*v*dT for agents in experiments

dT = 1/25; % Time step

[pos_a,length,times] = cut(file,agent); % Turns file into a position matrix, without NaN:s

myCircle = [pos_a(1,1,2000:size(pos_a,3)), pos_a(1,2,2000:size(pos_a,3))];
figure
x = myCircle(1,1,:);
y = myCircle(1,2,:);
x = x(:,:)';
y = y(:,:)';
plot(x,y,'c')
hold on
axis equal
%% Run to save results to file
clc;
% Sparar:
    % Filnamn/path (i en egen fil)
    % Kiralitet
    % NormArea/tid
    % Total tid
    % hastighet
    % l (size of area elements)
    % D_r

expName = 'circle_medium_1agent'; %Change name for each new set of data

file1 = ['results/Lab/' expName '.txt']; % Name of dataFile
file2 = ['results/Lab/' expName 'SourceFiles.txt']; % Name of file containing names of XML files
file3 = ['results/Lab/' expName 'indices.txt']

data = [result];
dlmwrite(file1,data,'-append');

fileID = fopen(file2,'a');
fprintf(fileID,'%-40s\n',file);
fclose(fileID);

dlmwrite(file3, indice, '-append')

%% Plot M, trajectory of agent

figure(3004)

x = myCircle(1,1,80:300);
y = myCircle(1,2,80:300);
x = x(:,:)';
y = y(:,:)';
plot(x,y,'c')
hold on
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

figure(2)
plot(v,'o')

figure(113)
semilogx(w,normArea,'o')
figure(114)
semilogx(w./v,normArea,'o')

%axis([0.01, 10, 0, 1.2])
%name=strcat('step', num2str(dT), '; ', 'time', num2str(numTimeSteps), '; ', 'simulations', num2str(numSimulations), '; ', obstacleType, '; ', 'R=', num2str(R), '; D_r=', num2str(D_r));
name = 'Experimental results';
title(name)
