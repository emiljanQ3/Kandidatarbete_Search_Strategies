%% PLot periodik results
clear,clf,clc
N_K = 15;

MarkSize = 40;
linewidth = 4;

tit_Font = 50;
ax_Font = 50;
%% experiment
clc
load('results/Final_results/2019425-1553_hm1agent_l216.mat')
clear 'length'

figure(118)


semilogx(abs(kir), normA, 'o')
title('rawdata')
axis([0.1 10 0 1.5])

figure(119)
hold on
[meanArea,binKir] = makeMean(kir,N_K,normA);
semilogx(binKir,meanArea,'k.','markersize',MarkSize)
title('Medelvärde over kiralitetbins')
axis([0.1 10 0 1.5])


%% sim per
load('results/Final_results/2019425-1345_hm_edge_t20.mat')

figure(119)
areaPerTime_max = numAgents*(4*v*measurmentTime*l/pi-l^2)/measurmentTime;
plot(w,areaPerTime/(areaPerTime_max),'r','Linewidth',linewidth)
axis([0.01, 10, 0, 1.2])
name=strcat('step', num2str(dT), '; ', 'time', num2str(numTimeSteps), '; ', 'simulations', num2str(numSimulations), '; ', obstacleType, '; ', 'R=', num2str(R), '; D_r=', num2str(D_r));
title(name)

%% sim beg
load('results/Final_results/2019425-1342_hm_infinite_t20.mat')

figure(119)
areaPerTime_max = numAgents*(4*v*measurmentTime*l/pi-l^2)/measurmentTime;
plot(w,areaPerTime/(areaPerTime_max),'b','Linewidth',linewidth)
axis([0.01, 10, 0, 1])
name=strcat('step', num2str(dT), '; ', 'time', num2str(numTimeSteps), '; ', 'simulations', num2str(numSimulations), '; ', obstacleType, '; ', 'R=', num2str(R), '; D_r=', num2str(D_r));
title(name)   
%% 

figure(119)
set(gca,'xscale','log')
axis([0.1 10 0 1.2])
title('Experimentell och simulerad data', 'Interpreter', 'latex')      %titla
set(gca, 'fontsize', tit_Font)

ylabel('Normerad area per tid', 'Interpreter', 'latex', 'fontsize', ax_Font)
xlabel('Kiralitet (rad/s)', 'Interpreter', 'latex', 'fontsize', ax_Font)




