%% PLot periodik results
clear,clf,clc
N_K = 15;

MarkSize = 30;
linewidth = 4;

tit_Font = 30;
ax_Font = 30;
num_Font = 20;

%% sim inf
load('results/Final_results/simulation/2019425-1342_hm_infinite_t20.mat')

figure(119)
areaPerTime_max = numAgents*(4*v*measurmentTime*l/pi-l^2)/measurmentTime;
P2 = plot(w,areaPerTime/(areaPerTime_max),'b','Linewidth',linewidth)
axis([0.01, 10, 0, 1])
 
%% experiment
clc
load('results/Final_results/experimental/2019425-1553_hm1agent_l216.mat')
clear 'length'

figure(118)


semilogx(abs(kir), normA, 'o')
title('rawdata')
axis([0.1 10 0 1.5])

figure(119)
hold on
[meanArea,binKir] = makeMean(kir,N_K,normA);
semilogx(binKir,meanArea,'k.','markersize',MarkSize+10)
P1 = plot(binKir,meanArea,'r.','markersize',MarkSize)
axis([0.1 10 0 1.5])

 
%% 
figure(119)
set(gca,'xscale','log')
axis([0 10 0 1.2])
set(gca, 'fontsize', num_Font)


xlabel('Kiralitet (rad/s)', 'Interpreter', 'latex', 'fontsize', ax_Font)
xticks([0.1 1 10])
xticklabels({'0.1','1','10'})


ylabel('S\"okeffektivitet','Position',[0.1 1.23], 'Interpreter','latex','fontsize', ax_Font)
set(get(gca,'ylabel'),'rotation',0)

legend('simulering','Location','northeast')
%legend([P2 P1],'simulering','experiment','Location','northeast')

set(gcf,'position',[100 100 900 800])
axis('square')
%%
%print(119,'Homogen_one','-djpeg')
%print(119,'Homogen_both','-djpeg')



