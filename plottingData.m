%% periodic
load('results/noEdge/oneAgent/homogeneous/2019218-92_hm02005_1.mat')

    figure(112)
    meanArea_max = meanAreaCovered(1);
    semilogx(w,meanAreaCovered/(meanArea_max),'*')
   % axis([0.01, 10, 0, 1.2])
    name=strcat('step', num2str(dT), '; ', 'time', num2str(numTimeSteps), '; ', 'simulations', num2str(numSimulations), '; ', obstacleType, '; ', 'R=', num2str(R), '; D_r=', num2str(D_r));
    title(name)
    hold on
    figure(113)
    plot(w,meanAreaCovered/(meanArea_max),'*')
    hold on
    
%% with edge
load('results/Edge/oneAgent/201937-1428_hm060167_1.mat')
if(edge) 
    figure(112)
    areaPerTime_max = numAgents*(4*v*measurmentTime*l/pi-l^2)/measurmentTime;
    semilogx(w,areaPerTime/(areaPerTime_max),'.')
    %axis([0.01, 10, 0, 1])
    name=strcat('step', num2str(dT), '; ', 'time', num2str(numTimeSteps), '; ', 'simulations', num2str(numSimulations), '; ', obstacleType, '; ', 'R=', num2str(R), '; D_r=', num2str(D_r));
    title(name)
    figure(113)
    plot(w,areaPerTime/(areaPerTime_max),'.')
else

    meanArea_max = numAgents*(4*v*measurmentTime*l/pi-l^2);
    semilogx(w,meanAreaCovered/(meanArea_max),'.')
    axis([0.01, 10, 0, 1.2])
    name=strcat('step', num2str(dT), '; ', 'time', num2str(numTimeSteps), '; ', 'simulations', num2str(numSimulations), '; ', obstacleType, '; ', 'R=', num2str(R), '; D_r=', num2str(D_r));
    title(name)
end

%%
experiment = 'hm1agent';
dataFile = ['results/Lab/' experiment '.txt'];
sourceFile = ['results/Lab/' experiment 'SourceFiles.txt'];
allData = dlmread(dataFile);

w = allData(:,1);
normArea = allData(:,2); 
v = allData(:,4);

figure(112)
semilogx(abs(w),normArea,'o')
figure(113)
plot(w,normArea,'o')

%axis([0.01, 10, 0, 1.2])
%name=strcat('step', num2str(dT), '; ', 'time', num2str(numTimeSteps), '; ', 'simulations', num2str(numSimulations), '; ', obstacleType, '; ', 'R=', num2str(R), '; D_r=', num2str(D_r));
name = 'Experimental results';
title(name)