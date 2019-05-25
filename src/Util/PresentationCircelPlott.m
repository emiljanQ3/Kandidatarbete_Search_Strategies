%% PLot cirkular results area over time
clf,clc
N_Kir = 15;

T       = 50;           % vid vilken tidpunkt plottar vi resultatet mot kiralitet
T2      = 50;           % Hur lång tid plottar vi arean 

index1   = floor(numAreaDP*T/measurmentTime);
index2   = floor(numAreaDP*T2/measurmentTime);

MarkSize = 30;
linewidth = 4;

tit_Font = 30;
ax_Font = 30;
num_Font = 20;
%%
%load('results/Final_results/simulation/2019417-1359_circle_R17_t360_l0156.mat')
%load('results/Final_results/simulation/2019425-943_circle_R058824_t180_l0156.mat')
load('results/Final_results/simulation/2019417-1057_circle_R1_t180_l0156.mat')
clear 'min'
clear 'max'
area_max = maxArea

str = 'b'

figure(5)
plot(1:index2, meanAreaCovered(1:index2,16)./maxArea,str,'Linewidth', linewidth)
hold on
figure(6)
plot(1:index2, meanAreaCovered(1:index2,50)./maxArea,str,'Linewidth', linewidth)
hold on
figure(7)
plot(1:index2, meanAreaCovered(1:index2,83)./maxArea,str,'Linewidth', linewidth)
hold on
%%
clc
%load('results/Final_results/experimental/2019424-1357_circle_large_1agent_t90_l216.mat')
%load('results/Final_results/experimental/2019426-111_circle_small2_1agent_t30_l216.mat')
load('results/Final_results/experimental/2019417-1959_circle_medium_1agent_t50_l216.mat')
clear 'length'

[meanArea,binKir, limit] = makeMean(kir,N_Kir,area);

str = 'r'

figure(5)
plot((1:size(meanArea,1))*T2/size(meanArea,1), meanArea(:,3)/maxArea,str,'Linewidth', linewidth)
hold on
figure(6)
hold on
plot((1:size(meanArea,1))*T2/size(meanArea,1),meanArea(:,9)/maxArea,str,'Linewidth', linewidth)
figure(7)
hold on
plot((1:size(meanArea,1))*T2/size(meanArea,1), meanArea(:,end)/maxArea,str,'Linewidth', linewidth)
%%
figure(8)
hold on

plot(binKir,meanArea(end,:)/maxArea,'k.','markersize',MarkSize+10)
P2 = plot(binKir,meanArea(end,:)/maxArea,'r.','markersize',MarkSize)
P1 = plot(w,meanAreaCovered(index1,:)/area_max,'b','Linewidth',linewidth)

axis(gca,[0 10 0 1]);
set(gca, 'fontsize', num_Font);
set(gca,'xscale','log')

xticks([0.1 1 10])
xticklabels({'0.1','1','10'})
xlabel('Kiralitet (rad/s)', 'Interpreter', 'latex', 'fontsize', ax_Font)

ylabel('Cirkelarea','Position',[0.1 1.01], 'Interpreter','latex','fontsize', ax_Font)
set(get(gca,'ylabel'),'rotation',0)
set(gcf,'position',[100 100 900 800])
axis(gca,'square')
legend([P1 P2],'simulering','experiment','Location','northwest')
%%
clc

allAxes = findall(0, 'Type', 'axes');
allFigures = findall(0, 'Type', 'figure');

axis(allAxes,[0 50 0 1]);
set(allAxes, 'fontsize', num_Font);
for i = 5:7
    figure(i)
    xlabel('Tid (s)', 'Interpreter', 'latex', 'fontsize', ax_Font)
  
    ylabel('Cirkelarea','Position',[0.1 1.01], 'Interpreter','latex','fontsize', ax_Font)
    set(get(gca,'ylabel'),'rotation',0)
    
    legend('simulering','experiment','Location','northwest')
end
set(allFigures,'position',[100 100 900 800])
axis(allAxes,'square')

%%
%print(5,'Circle_both_low','-djpeg')
%print(5,'Circle_one__low','-djpeg')

%print(6,'Circle_both_med','-djpeg')
%print(6,'Circle_one__med','-djpeg')

%print(7,'Circle_both_high','-djpeg')
%print(7,'Circle_one__high','-djpeg')

%print(8,'Circle_crossSection_medium','-djpeg')

