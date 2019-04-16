%% Simulated results

load('results/Final_simulation/2019416-144_circle_R1_t50_l0156.mat')
    T = 50;         % vid vilken tidpunkt plottar vi resultatet
    index = floor(numAreaDP*T/measurmentTime);

    figure(111)
    hold on
    c = jet(length(w));
    %c(:,2) = 0;
    for i = 1:5:size(meanAreaCovered,2)
        plot(1:numAreaDP, meanAreaCovered(:,i)./maxArea, 'color', c(i,:))
    end


    figure(112)
    hold on
    for i = 1:size(meanAreaCovered,2)
        plot(w(i), meanAreaCovered(index,i)/(pi*R^2),'o','color',c(i,:))
    end
    set(gca,'xscale','log')
    
    
    title('')
    name=strcat('step', num2str(dT), '; ', 'time', num2str(numTimeSteps), '; ', 'simulations', num2str(numSimulations), '; ', obstacleType, '; ', 'R=', num2str(R), '; D_r=', num2str(D_r));
    title(name)

%% Experimental results

experiment = 'hm1agent';
dataFile = ['results/Lab/' experiment '.txt'];
sourceFile = ['results/Lab/' experiment 'SourceFiles.txt'];
allData = dlmread(dataFile);

w = allData(:,1);
normArea = allData(:,2); 
v = allData(:,4);

figure(112)
semilogx(abs(w),normArea,'bo')
figure(113)
plot(w,normArea,'o')

%axis([0.01, 10, 0, 1.2])
%name=strcat('step', num2str(dT), '; ', 'time', num2str(numTimeSteps), '; ', 'simulations', num2str(numSimulations), '; ', obstacleType, '; ', 'R=', num2str(R), '; D_r=', num2str(D_r));
name = 'Experimental results';
title(name)

%%

figure(112)

axis([0.02 10 0 1.3])
title('Av agenter uppsökt area normerat mot maximal upptäckt area', 'Interpreter', 'latex')      %titla
set(gca, 'fontsize', 28)

ylabel('Effektivitet', 'Interpreter', 'latex', 'fontsize', 35)
xlabel('Kiralitet (rad/s)', 'Interpreter', 'latex', 'fontsize', 35)

legend('Periodisk','Begränsad','Exprimentel')



